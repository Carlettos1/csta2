use proc_macro2::{Span, TokenStream};
use quote::{quote, quote_spanned};
use rand::distr::weighted;
use syn::spanned::Spanned;
use syn::*;

#[proc_macro_derive(Randomizable, attributes(csta))]
pub fn derive_randomizable(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = input.ident;

    let generics = add_trait_bounds(input.generics);
    let (impl_generics, ty_generics, where_clause) = generics.split_for_impl();
    match input.data {
        Data::Struct(data) => match data.fields {
            Fields::Named(fields) => {
                let (let_quotes, field_quotes) = parse_fields_named(&fields);
                quote! {
                    impl #impl_generics csta::Randomizable for #name #ty_generics #where_clause {
                        #[allow(unused)]
                        fn sample<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
                            #( #let_quotes; )*
                            Self {
                                #( #field_quotes, )*
                            }
                        }
                    }
                }
            }
            Fields::Unnamed(fields) => {
                let random_fields = parse_fields_unnamed(&fields);
                quote! {
                    impl #impl_generics csta::Randomizable for #name #ty_generics #where_clause {
                        fn sample<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
                            Self(
                                #( #random_fields, )*
                            )
                        }
                    }
                }
            }
            Fields::Unit => {
                quote! {
                    impl #impl_generics csta::Randomizable for #name #ty_generics #where_clause {
                        fn sample<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
                            Self
                        }
                    }
                }
            }
        },
        // todo: add weighted probabilities
        Data::Enum(data) => {
            // see if they have weighted probabilities
            if data.variants.iter().any(enum_has_attribute) {
                // at least one have weighted probabilities.
                // if one have it, then all must have it as well.
                assert!(
                    data.variants.iter().all(enum_has_attribute),
                    "If one variant has the weight attribute, all should.\nHint: add #[csta(weight = 0.1)] to ALL variants"
                );
                // I need, total_prob = SUM(weight)
                // if total_prob == 0.0, return default or first.
                // prob = weight / total_prob
                // r = rng()
                // if r < prob1 {1} else if r - prob1 < prob2 {2}

                let len_variants = data.variants.len();

                let probabilities = data.variants.iter().map(|variant| {
                    let enum_attributes = get_parsed_enum_attributes(variant);
                    #[allow(clippy::infallible_destructuring_match)]
                    let weight = match &enum_attributes[0] {
                        CstaEnumAttributes::Weighted(float) => float,
                    };
                    
                    quote_spanned! {variant.span()=>
                        #weight
                    }
                });

                let builders = data.variants.iter().map(|variant| {
                    let iden = &variant.ident;
                    match &variant.fields {
                        Fields::Named(fields) => {
                            let (let_quotes, field_quotes) = parse_fields_named(fields);
                            quote_spanned! {variant.span()=>
                                {
                                    #( #let_quotes; )*
                                    #name::#iden { #( #field_quotes, )* }
                                }
                            }
                        }
                        Fields::Unnamed(fields) => {
                            let random_fields = parse_fields_unnamed(fields);
                            quote_spanned! {variant.span()=>
                                #name::#iden( #( #random_fields, )* )
                            }
                        }
                        Fields::Unit => {
                            quote_spanned! {variant.span()=>
                                #name::#iden
                            }
                        }
                    }
                }).collect::<Vec<_>>();

                let default = &builders[0];
                let probabilities: Vec<_> = probabilities.into_iter().zip(data.variants.iter()).scan(quote!(0.0_f64), |state, (prob, variant)| {
                    let tmp = quote_spanned! {variant.span()=>
                        #state + #prob
                    };
                    *state = tmp;
                    Some(state.clone())
                }).collect();

                let prob_sum = probabilities.last().unwrap();

                let if_builder_chain = probabilities.iter().zip(builders.iter()).map(|(prob, builder)| {
                    quote_spanned! {prob.span()=>
                        if r < #prob {
                            return #builder;
                        }
                    }
                });

                quote! {
                    impl #impl_generics csta::Randomizable for #name #ty_generics #where_clause {
                        #[allow(unused)]
                        fn sample<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
                            let total_probability = #prob_sum;
                            if total_probability == 0.0 {
                                return #default;
                            }

                            let mut r: f64 = rng.random::<f64>() * total_probability;
                            #( #if_builder_chain )*

                            #default
                        }
                    }
                }
            } else {
                // if no one have weighted, just use the N-dice approach.
                let num = data.variants.len();
                let random_variants = data.variants.iter().enumerate().map(|(i, variant)| {
                    let index = Index::from(i);
                    let iden = &variant.ident;

                    match &variant.fields {
                        Fields::Named(fields) => {
                            let (let_quotes, field_quotes) = parse_fields_named(fields);
                            quote_spanned! {variant.span()=>
                                #index => {
                                    #( #let_quotes; )*
                                    #name::#iden { #( #field_quotes, )* }
                                }
                            }
                        }
                        Fields::Unnamed(fields) => {
                            let random_fields = parse_fields_unnamed(fields);
                            quote_spanned! {variant.span()=>
                                #index => #name::#iden( #( #random_fields, )* )
                            }
                        }
                        Fields::Unit => {
                            quote_spanned! {variant.span()=>
                                #index => #name::#iden
                            }
                        }
                    }
                });
                quote! {
                    impl #impl_generics csta::Randomizable for #name #ty_generics #where_clause {
                        #[allow(unused)]
                        fn sample<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
                            let num = rng.random_range(0..#num);
                            match num {
                                #( #random_variants, )*
                                _ => unreachable!("Number not in range of enum"),
                            }
                        }
                    }
                }
            }
        }
        Data::Union(_) => unimplemented!(),
    }
    .into()
}

fn add_trait_bounds(mut generics: Generics) -> Generics {
    for param in &mut generics.params {
        if let GenericParam::Type(ref mut type_param) = *param {
            type_param.bounds.push(parse_quote!(csta::Randomizable));
        }
    }
    generics
}

enum CstaEnumAttributes {
    Weighted(LitFloat),
}

fn enum_has_attribute(variant: &Variant) -> bool {
    let mut csta_attributes = Vec::new();
    parse_enum_attributes(&variant.attrs, &mut csta_attributes);
    !csta_attributes.is_empty()
}

fn get_parsed_enum_attributes(variant: &Variant) -> Vec<CstaEnumAttributes> {
    let mut csta_attributes = Vec::new();
    parse_enum_attributes(&variant.attrs, &mut csta_attributes);
    csta_attributes
}

fn parse_enum_attributes(
    attributes: &Vec<Attribute>,
    csta_attributes: &mut Vec<CstaEnumAttributes>,
) {
    for attr in attributes {
        if attr.path().is_ident("csta") {
            attr.parse_nested_meta(|meta| {
                if meta.path.is_ident("weight") {
                    if let Ok(value) = meta.value() {
                        let expr: Expr = value.parse()?;
                        if let Expr::Lit(lit) = expr {
                            if let Lit::Float(float) = lit.lit {
                                csta_attributes.push(CstaEnumAttributes::Weighted(float));
                            } else {
                                return Err(Error::new(attr.span(), "Expected a float number"));
                            }
                        } else {
                            return Err(Error::new(attr.span(), "Expected a float number"));
                        }
                    } else {
                        return Err(Error::new(attr.span(), "Expected a float number"));
                    }
                }
                Ok(())
            })
            .expect("Failed to parse attribute");
        }
    }
}

fn parse_attribute(attr: &Attribute, csta_attribute: &mut CstaAttributes) {
    if attr.path().is_ident("csta") {
        attr.parse_nested_meta(|meta| {
            if meta.path.is_ident("range") {
                let content;
                parenthesized!(content in meta.input);
                let range: Expr = content.parse()?;
                if let Expr::Range(range) = range {
                    // Check that the range has start and end
                    if range.start.is_none() || range.end.is_none() {
                        return Err(Error::new(
                            range.span(),
                            "Expected range with start and end (either a..b or a..=b)",
                        ));
                    }
                    *csta_attribute = CstaAttributes::Range(range);
                } else {
                    return Err(Error::new(
                        range.span(),
                        "Expected range (either a..b or a..=b)",
                    ));
                }
            }
            if meta.path.is_ident("len") {
                let content;
                parenthesized!(content in meta.input);
                let expr: Expr = content.parse()?;
                *csta_attribute = CstaAttributes::Len(expr);
            }
            if meta.path.is_ident("after") {
                let content;
                parenthesized!(content in meta.input);
                let expr: Expr = content.parse()?;
                *csta_attribute = CstaAttributes::After(expr);
            }
            if meta.path.is_ident("default") {
                if let Ok(value) = meta.value() {
                    let iden: TokenStream = value.parse()?;
                    *csta_attribute = CstaAttributes::DefaultWith(iden);
                } else {
                    *csta_attribute = CstaAttributes::Default;
                }
            }
            if meta.path.is_ident("mul") {
                let value = meta.value()?;
                csta_attribute.add_mul(Mul(value.parse()?));
            }
            if meta.path.is_ident("div") {
                let value = meta.value()?;
                csta_attribute.add_div(Div(value.parse()?));
            }
            if meta.path.is_ident("add") {
                let value = meta.value()?;
                csta_attribute.add_add(Add(value.parse()?));
            }
            if meta.path.is_ident("sub") {
                let value = meta.value()?;
                csta_attribute.add_sub(Sub(value.parse()?));
            }
            Ok(())
        })
        .expect("Failed to parse attribute");
    }
}

enum CstaAttributes {
    UseRandomizable,
    Range(ExprRange),
    Len(Expr), // used in Vec<T>
    // TODO: Probability(TokenStream), // used in Option<T>
    After(Expr), // for manipulations after being created via randomizable
    Default,
    DefaultWith(TokenStream),
    Operation(Option<Mul>, Option<Div>, Option<Add>, Option<Sub>),
}

impl CstaAttributes {
    pub fn set_op(&mut self) {
        if !matches!(self, CstaAttributes::Operation(_, _, _, _)) {
            *self = CstaAttributes::Operation(None, None, None, None);
        }
    }

    pub fn add_mul(&mut self, value: Mul) {
        self.set_op();
        if let CstaAttributes::Operation(mul, _, _, _) = self {
            *mul = Some(value);
        }
    }

    pub fn add_div(&mut self, value: Div) {
        self.set_op();
        if let CstaAttributes::Operation(_, div, _, _) = self {
            *div = Some(value);
        }
    }

    pub fn add_add(&mut self, value: Add) {
        self.set_op();
        if let CstaAttributes::Operation(_, _, add, _) = self {
            *add = Some(value);
        }
    }

    pub fn add_sub(&mut self, value: Sub) {
        self.set_op();
        if let CstaAttributes::Operation(_, _, _, sub) = self {
            *sub = Some(value);
        }
    }
}

struct Mul(TokenStream);
struct Div(TokenStream);
struct Add(TokenStream);
struct Sub(TokenStream);

// the different thing between named and unnamed is that named fields should be able to be used on other attributes
// like #[csta(len(w*h))], with w, h being other fields.
// in unnamed fields is imposible to do that, so its a different bussiness
/// returns (let_quotes, fields_quotes), in that order
fn parse_fields_named(fields: &FieldsNamed) -> (Vec<TokenStream>, Vec<TokenStream>) {
    let mut early_let_quotes = Vec::new();
    let mut later_let_quotes = Vec::new();
    let mut last_let_quotes = Vec::new();
    let mut fields_quotes = Vec::new();

    for field in &fields.named {
        let mut attribute = CstaAttributes::UseRandomizable; // w/o attributes, use randomizable as default
        field
            .attrs
            .iter()
            .for_each(|attr| parse_attribute(attr, &mut attribute));
        let ident = &field.ident;
        let field_type = &field.ty;
        let value = apply_attributes(field_type, field.span(), &attribute);
        match attribute {
            CstaAttributes::Default => {
                // Default::default() will get the earlier priority.
                early_let_quotes.push(quote_spanned! {field.span()=>
                    let #ident = #value
                });
            }
            CstaAttributes::DefaultWith(_) => {
                // These will get the second priority, so that they can use default fields
                later_let_quotes.push(quote_spanned! {field.span()=>
                    let #ident = #value
                });
            }
            CstaAttributes::After(expr) => {
                // after only works on named for now.
                // it creates a let = T::sample, and then a let = expr;
                later_let_quotes.push(quote_spanned! {field.span()=>
                    let #ident = <#field_type as ::csta::Randomizable>::sample(rng)
                });
                last_let_quotes.push(quote_spanned! {field.span()=>
                    let #ident = #expr
                });
            }
            _ => {
                // These are last prio, maybe they are in order so their prio is in written order
                last_let_quotes.push(quote_spanned! {fields.span()=>
                    let #ident = #value
                });
            }
        }
        // because everything is a let w = #value, Self { w } is used.
        fields_quotes.push(quote_spanned! {fields.span()=>
            #ident
        });
    }
    // now we merge let quotes and return everything in correct order.
    early_let_quotes.append(&mut later_let_quotes);
    early_let_quotes.append(&mut last_let_quotes);
    (early_let_quotes, fields_quotes)
}

fn parse_fields_unnamed(fields: &FieldsUnnamed) -> impl Iterator<Item = TokenStream> + '_ {
    fields.unnamed.iter().map(|field| {
        let mut modifier = CstaAttributes::UseRandomizable;
        field
            .attrs
            .iter()
            .for_each(|attr| parse_attribute(attr, &mut modifier));

        let field_type = &field.ty;
        apply_attributes(field_type, field.span(), &modifier)
    })
}

fn apply_attributes(field_type: &Type, span: Span, modifier: &CstaAttributes) -> TokenStream {
    match modifier {
        CstaAttributes::UseRandomizable => quote_spanned! {span=>
            <#field_type as ::csta::Randomizable>::sample(rng)
        },
        CstaAttributes::Range(range) => quote_spanned! {span=>
            rng.random_range(#range)
        },
        CstaAttributes::Default => quote_spanned! {span=>
            Default::default()
        },
        CstaAttributes::DefaultWith(iden) => quote_spanned! {span=>
            #iden
        },
        CstaAttributes::Len(len) => {
            // if is_vec(field_type) {
            //     let generics = extract_vec_inner(field_type);
            //     if let Some(inner_type) = generics {
            //         apply_modifier(inner_type, span, modifier)
            //     } else {
            //         panic!("Vec needs to have generics (Vec<T>)");
            //     }
            // } else {
            //     quote_spanned! (span=>)
            // }
            let generics = extract_vec_inner(field_type);
            if let Some(inner_type) = generics {
                quote_spanned! {span=>
                    (0..#len).map(|_| <#inner_type as ::csta::Randomizable>::sample(rng)).collect()
                }
            } else {
                quote_spanned! (span=>)
            }
        }
        CstaAttributes::After(expr) => quote_spanned! {span=>
            #expr
        },
        CstaAttributes::Operation(mul, div, add, sub) => {
            let mut field = quote_spanned! {span=>
                #field_type::sample(rng)
            };
            if let Some(Mul(mul)) = mul {
                field.extend(quote_spanned! {span=>
                    * #mul
                });
            }
            if let Some(Div(div)) = div {
                field.extend(quote_spanned! {span=>
                    / #div
                });
            }
            if let Some(Add(add)) = add {
                field.extend(quote_spanned! {span=>
                    + #add
                });
            }
            if let Some(Sub(sub)) = sub {
                field.extend(quote_spanned! {span=>
                    - #sub
                });
            }
            field
        }
    }
}

fn extract_vec_inner(ty: &Type) -> Option<&Type> {
    if let Type::Path(type_path) = ty
        && let Some(last_segment) = type_path.path.segments.last()
        && last_segment.ident == "Vec"
        && let PathArguments::AngleBracketed(ref generic_args) = last_segment.arguments
        && let Some(GenericArgument::Type(inner_ty)) = generic_args.args.first()
    {
        return Some(inner_ty);
    }
    None
}

fn is_vec(ty: &Type) -> bool {
    if let Type::Path(type_path) = ty
        && let Some(last_segment) = type_path.path.segments.last()
        && last_segment.ident == "Vec"
    {
        true
    } else {
        false
    }
}
