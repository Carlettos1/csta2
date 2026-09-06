#![doc = include_str!("../README.md")]
use proc_macro2::{Span, TokenStream};
use quote::{quote, quote_spanned};
use syn::{spanned::Spanned, *};

#[proc_macro_derive(Randomizable, attributes(csta))]
pub fn derive_randomizable(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    match expand(parse_macro_input!(input as DeriveInput)) {
        Ok(t) => t,
        Err(e) => e.to_compile_error(),
    }
    .into()
}
fn expand(input: DeriveInput) -> Result<TokenStream> {
    let name = input.ident;
    let mut generics = input.generics;
    let body = match input.data {
        Data::Struct(s) => fields(&s.fields, quote!(Self), &mut generics)?,
        Data::Union(u) => {
            return Err(Error::new(
                u.union_token.span,
                "Randomizable does not support unions",
            ));
        }
        Data::Enum(e) => {
            if e.variants.is_empty() {
                return Err(Error::new(name.span(), "cannot sample an empty enum"));
            }
            let mut weights = Vec::new();
            let mut builders = Vec::new();
            for v in &e.variants {
                weights.push(weight(v)?);
                let id = &v.ident;
                builders.push(fields(&v.fields, quote!(Self::#id), &mut generics)?);
            }
            if weights.iter().any(Option::is_some) {
                if weights.iter().any(Option::is_none) {
                    return Err(Error::new(
                        name.span(),
                        "all variants must specify a weight",
                    ));
                }
                let raw: Vec<_> = weights.into_iter().map(Option::unwrap).collect();
                let max = raw.iter().copied().fold(0.0_f64, f64::max);
                if max == 0.0 {
                    return Err(Error::new(
                        name.span(),
                        "at least one weight must be positive",
                    ));
                }
                // Scale first to avoid overflow of sums of valid finite weights.
                let total: f64 = raw.iter().map(|x| x / max).sum();
                let mut cumulative = 0.0;
                let mut branches = Vec::new();
                let mut last = None;
                for (w, builder) in raw.iter().zip(&builders) {
                    if *w == 0.0 {
                        continue;
                    }
                    cumulative += (*w / max) / total;
                    let threshold = cumulative;
                    branches.push(quote!(if draw < #threshold {return #builder;}));
                    last = Some(builder);
                }
                let fallback = last.unwrap();
                quote!({let draw=rand::RngExt::random::<f64>(rng);#(#branches)* #fallback})
            } else {
                let n = builders.len();
                let indexes = 0..n;
                quote!({match rand::RngExt::random_range(rng,0..#n) {#(#indexes=>#builders,)* _=>unreachable!("range checked")}})
            }
        }
    };
    let (ig, tg, wc) = generics.split_for_impl();
    Ok(quote!(impl #ig csta::Randomizable for #name #tg #wc {
        #[allow(unused_variables,unused_parens)]
        fn sample<R:rand::RngExt+?Sized>(rng:&mut R)->Self {#body}
    }))
}
fn weight(v: &Variant) -> Result<Option<f64>> {
    let mut result = None;
    for attr in &v.attrs {
        if attr.path().is_ident("csta") {
            attr.parse_nested_meta(|m| {
                if !m.path.is_ident("weight") {
                    return Err(m.error("unknown variant attribute"));
                }
                if result.is_some() {
                    return Err(m.error("duplicate weight"));
                }
                let expr: Expr = m.value()?.parse()?;
                let n = literal_number(&expr).ok_or_else(|| {
                    Error::new(
                        expr.span(),
                        "weight must be a finite nonnegative numeric literal",
                    )
                })?;
                if !n.is_finite() || n < 0.0 {
                    return Err(Error::new(
                        expr.span(),
                        "weight must be finite and nonnegative",
                    ));
                }
                result = Some(n);
                Ok(())
            })?;
        }
    }
    Ok(result)
}
fn literal_number(e: &Expr) -> Option<f64> {
    match e {
        Expr::Lit(e) => match &e.lit {
            Lit::Float(x) => x.base10_parse().ok(),
            Lit::Int(x) => x.base10_parse().ok(),
            _ => None,
        },
        Expr::Unary(e) if matches!(e.op, UnOp::Neg(_)) => literal_number(&e.expr).map(|x| -x),
        _ => None,
    }
}
#[derive(Default)]
struct Attr {
    base: Option<Base>,
    ops: Vec<(String, Expr)>,
}
enum Base {
    Range(ExprRange),
    Len(Expr),
    Default(Option<Expr>),
    After(Expr),
}
fn attrs(f: &Field) -> Result<Attr> {
    let mut a = Attr::default();
    for attr in &f.attrs {
        if attr.path().is_ident("csta") {
            attr.parse_nested_meta(|m| {
                let id = m
                    .path
                    .get_ident()
                    .ok_or_else(|| m.error("unknown field attribute"))?
                    .to_string();
                if ["mul", "div", "add", "sub"].contains(&id.as_str()) {
                    if a.base.is_some() || a.ops.iter().any(|(i, _)| *i == id) {
                        return Err(m.error("conflicting or duplicate field attribute"));
                    }
                    let e: Expr = m.value()?.parse()?;
                    if id == "div" && literal_number(&e) == Some(0.0) {
                        return Err(m.error("division by zero"));
                    }
                    a.ops.push((id, e));
                    return Ok(());
                }
                if a.base.is_some() || !a.ops.is_empty() {
                    return Err(m.error("conflicting field initialization attributes"));
                }
                a.base = Some(match id.as_str() {
                    "default" => Base::Default(if m.input.peek(Token![=]) {
                        Some(m.value()?.parse()?)
                    } else {
                        None
                    }),
                    "range" | "len" | "after" => {
                        let content;
                        parenthesized!(content in m.input);
                        let e: Expr = content.parse()?;
                        match id.as_str() {
                            "range" => {
                                let Expr::Range(r) = e else {
                                    return Err(m.error("expected bounded range"));
                                };
                                let (Some(lo), Some(hi)) = (&r.start, &r.end) else {
                                    return Err(m.error("range needs both endpoints"));
                                };
                                if let (Some(l), Some(h)) = (literal_number(lo), literal_number(hi))
                                    && (!l.is_finite()
                                        || !h.is_finite()
                                        || l > h
                                        || (l == h && matches!(r.limits, RangeLimits::HalfOpen(_))))
                                {
                                    return Err(m.error("range is empty, reversed or nonfinite"));
                                }
                                Base::Range(r)
                            }
                            "len" => {
                                if inner_vec(&f.ty).is_none() {
                                    return Err(m.error("len requires Vec<T>"));
                                }
                                if literal_number(&e)
                                    .is_some_and(|n| n < 0.0 || n.fract() != 0.0 || !n.is_finite())
                                {
                                    return Err(m.error("length must be a nonnegative integer"));
                                }
                                Base::Len(e)
                            }
                            _ => Base::After(e),
                        }
                    }
                    _ => return Err(m.error("unknown field attribute")),
                });
                Ok(())
            })?;
        }
    }
    Ok(a)
}
fn inner_vec(t: &Type) -> Option<&Type> {
    if let Type::Path(p) = t
        && let Some(s) = p.path.segments.last()
        && s.ident == "Vec"
        && let PathArguments::AngleBracketed(a) = &s.arguments
        && let Some(GenericArgument::Type(t)) = a.args.first()
    {
        Some(t)
    } else {
        None
    }
}
fn sample_type(t: &Type, g: &mut Generics) -> TokenStream {
    g.make_where_clause()
        .predicates
        .push(parse_quote!(#t:csta::Randomizable));
    quote!(<#t as csta::Randomizable>::sample(rng))
}
fn fields(fs: &Fields, ctor: TokenStream, g: &mut Generics) -> Result<TokenStream> {
    let names: Vec<_> = fs
        .iter()
        .enumerate()
        .map(|(i, f)| {
            f.ident
                .clone()
                .unwrap_or_else(|| Ident::new(&format!("__csta_field_{i}"), Span::call_site()))
        })
        .collect();
    let mut early = Vec::new();
    let mut nodes = Vec::new();
    for (f, id) in fs.iter().zip(&names) {
        let t = &f.ty;
        let a = attrs(f)?;
        let mut after = false;
        let expr = match a.base {
            Some(Base::Range(r)) => quote_spanned!(f.span()=>rand::RngExt::random_range(rng,#r)),
            Some(Base::Len(e)) => {
                let sample = sample_type(inner_vec(t).unwrap(), g);
                quote!({let __csta_len=usize::try_from(#e).expect("vector length must fit usize");(0..__csta_len).map(|_|#sample).collect::<#t>()})
            }
            Some(Base::Default(None)) => {
                g.make_where_clause()
                    .predicates
                    .push(parse_quote!(#t:Default));
                quote!(Default::default())
            }
            Some(Base::Default(Some(e))) => quote!(#e),
            Some(Base::After(e)) => {
                let sample = sample_type(t, g);
                early.push(quote!(let #id:#t=#sample;));
                after = true;
                quote!(#e)
            }
            None => {
                let mut ex = sample_type(t, g);
                for (op, e) in a.ops {
                    ex = match op.as_str() {
                        "mul" => quote!((#ex)*(#e)),
                        "div" => quote!((#ex)/(#e)),
                        "add" => quote!((#ex)+(#e)),
                        _ => quote!((#ex)-(#e)),
                    };
                }
                ex
            }
        };
        struct Paths(Vec<String>);
        impl<'a> syn::visit::Visit<'a> for Paths {
            fn visit_expr_path(&mut self, p: &'a ExprPath) {
                if p.qself.is_none() && p.path.segments.len() == 1 {
                    self.0.push(p.path.segments[0].ident.to_string());
                }
            }
        }
        let parsed: Expr = syn::parse2(expr.clone())?;
        let mut paths = Paths(Vec::new());
        syn::visit::Visit::visit_expr(&mut paths, &parsed);
        let deps: Vec<usize> = names
            .iter()
            .enumerate()
            .filter(|(_, n)| paths.0.contains(&n.to_string()) && !(after && *n == id))
            .map(|(i, _)| i)
            .collect();
        nodes.push((deps, quote_spanned!(f.span()=>let #id:#t=#expr;)));
    }
    let mut done = vec![false; nodes.len()];
    let mut ordered = Vec::new();
    while ordered.len() < nodes.len() {
        let Some(i) = (0..nodes.len()).find(|i| !done[*i] && nodes[*i].0.iter().all(|d| done[*d]))
        else {
            return Err(Error::new(
                fs.span(),
                "cyclic field initialization dependencies",
            ));
        };
        done[i] = true;
        ordered.push(nodes[i].1.clone());
    }
    let value = match fs {
        Fields::Named(_) => quote!(#ctor {#(#names),*}),
        Fields::Unnamed(_) => quote!(#ctor (#(#names),*)),
        Fields::Unit => ctor,
    };
    Ok(quote!({#(#early)* #(#ordered)* #value}))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn diagnostics() {
        for input in [
            "enum E {}",
            "union U { x: f64 }",
            "enum E {#[csta(weight=0.0)] A}",
            "enum E {#[csta(weight=-1.0)] A}",
            "enum E {#[csta(weight=1e999)] A}",
            "enum E {#[csta(weight=1.0)] A,B}",
            "struct S {#[csta(foo)] x:f64}",
            "struct S {#[csta(range(1..1))] x:f64}",
            "struct S {#[csta(range(2..1))] x:f64}",
            "struct S {#[csta(range(..1))] x:f64}",
            "struct S {#[csta(len(2))] x:f64}",
            "struct S {#[csta(div=0)] x:f64}",
        ] {
            let ast = syn::parse_str(input).unwrap();
            assert!(expand(ast).is_err(), "{input}");
        }
        for input in [
            "struct S;",
            "enum E {#[csta(weight=0.0)] A,#[csta(weight=3.0)] B}",
            "struct S<T> {x:T}",
            "struct S {#[csta(len(0))] x:Vec<f64>}",
        ] {
            assert!(expand(syn::parse_str(input).unwrap()).is_ok());
        }
    }
}
