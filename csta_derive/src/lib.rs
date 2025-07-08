use proc_macro2::TokenStream;
use quote::{quote, quote_spanned};
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
                let random_fields = parse_fields_named(&fields);
                quote! {
                    impl #impl_generics csta::Randomizable for #name #ty_generics #where_clause {
                        fn sample<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
                            Self {
                                #( #random_fields, )*
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
        Data::Enum(data) => {
            let num = data.variants.len();
            let random_variants = data.variants.iter().enumerate().map(|(i, variant)| {
                let index = Index::from(i);
                let iden = &variant.ident;
                match &variant.fields {
                    Fields::Named(fields) => {
                        let random_fields = parse_fields_named(fields);
                        quote_spanned! {variant.span()=>
                            #index => #name::#iden { #( #random_fields, )* }
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
                    fn sample<R: rand::Rng + ?Sized>(rng: &mut R) -> Self {
                        let num = rng.gen_range(0..#num);
                        match num {
                            #( #random_variants, )*
                            _ => panic!("Number not in range of enum"),
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

fn parse_attribute(attr: &Attribute, random_type: &mut RandomField) {
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
                    *random_type = RandomField::Range(range);
                } else {
                    return Err(Error::new(
                        range.span(),
                        "Expected range (either a..b or a..=b)",
                    ));
                }
            }
            if meta.path.is_ident("default") {
                if let Ok(value) = meta.value() {
                    let iden: TokenStream = value.parse()?;
                    *random_type = RandomField::DefaultWith(iden);
                } else {
                    *random_type = RandomField::Default;
                }
            }
            if meta.path.is_ident("mul") {
                let value = meta.value()?;
                random_type.add_mul(Mul(value.parse()?));
            }
            if meta.path.is_ident("div") {
                let value = meta.value()?;
                random_type.add_div(Div(value.parse()?));
            }
            if meta.path.is_ident("add") {
                let value = meta.value()?;
                random_type.add_add(Add(value.parse()?));
            }
            if meta.path.is_ident("sub") {
                let value = meta.value()?;
                random_type.add_sub(Sub(value.parse()?));
            }
            Ok(())
        })
        .expect("Failed to parse attribute");
    }
}

enum RandomField {
    UseRandomizable,
    Range(ExprRange),
    Default,
    DefaultWith(TokenStream),
    Operation(Option<Mul>, Option<Div>, Option<Add>, Option<Sub>),
}

impl RandomField {
    pub fn set_op(&mut self) {
        if !matches!(self, RandomField::Operation(_, _, _, _)) {
            *self = RandomField::Operation(None, None, None, None);
        }
    }

    pub fn add_mul(&mut self, value: Mul) {
        self.set_op();
        if let RandomField::Operation(mul, _, _, _) = self {
            *mul = Some(value);
        }
    }

    pub fn add_div(&mut self, value: Div) {
        self.set_op();
        if let RandomField::Operation(_, div, _, _) = self {
            *div = Some(value);
        }
    }

    pub fn add_add(&mut self, value: Add) {
        self.set_op();
        if let RandomField::Operation(_, _, add, _) = self {
            *add = Some(value);
        }
    }

    pub fn add_sub(&mut self, value: Sub) {
        self.set_op();
        if let RandomField::Operation(_, _, _, sub) = self {
            *sub = Some(value);
        }
    }
}

struct Mul(TokenStream);
struct Div(TokenStream);
struct Add(TokenStream);
struct Sub(TokenStream);

fn parse_fields_named(fields: &FieldsNamed) -> impl Iterator<Item = TokenStream> + '_ {
    fields.named.iter().map(|field| {
        let mut random_type = RandomField::UseRandomizable;
        field
            .attrs
            .iter()
            .for_each(|attr| parse_attribute(attr, &mut random_type));
        let ident = &field.ident;
        let field_type = &field.ty;
        match random_type {
            RandomField::UseRandomizable => quote_spanned! {field.span()=>
                #ident: #field_type::sample(rng)
            },
            RandomField::Range(range) => quote_spanned! {field.span()=>
                #ident: rng.gen_range(#range)
            },
            RandomField::Default => quote_spanned! {field.span()=>
                #ident: Default::default()
            },
            RandomField::Operation(mul, div, add, sub) => {
                let mut field = quote_spanned! {field.span()=>
                    #ident: #field_type::sample(rng)
                };
                if let Some(Mul(mul)) = mul {
                    field.extend(quote_spanned! {field.span()=>
                        * #mul
                    });
                }
                if let Some(Div(div)) = div {
                    field.extend(quote_spanned! {field.span()=>
                        / #div
                    });
                }
                if let Some(Add(add)) = add {
                    field.extend(quote_spanned! {field.span()=>
                        + #add
                    });
                }
                if let Some(Sub(sub)) = sub {
                    field.extend(quote_spanned! {field.span()=>
                        - #sub
                    });
                }
                field
            }
            RandomField::DefaultWith(iden) => quote_spanned! {field.span()=>
                #ident: #iden
            },
        }
    })
}

fn parse_fields_unnamed(fields: &FieldsUnnamed) -> impl Iterator<Item = TokenStream> + '_ {
    fields.unnamed.iter().map(|field| {
        let mut random_type = RandomField::UseRandomizable;
        field
            .attrs
            .iter()
            .for_each(|attr| parse_attribute(attr, &mut random_type));

        let field_type = &field.ty;
        match random_type {
            RandomField::UseRandomizable => quote_spanned! {field.span()=>
                #field_type::sample(rng)
            },
            RandomField::Range(range) => quote_spanned! {field.span()=>
                rng.gen_range(#range)
            },
            RandomField::Default => quote_spanned! {field.span()=>
                Default::default()
            },
            RandomField::Operation(mul, div, add, sub) => {
                let mut field = quote_spanned! {field.span()=>
                    #field_type::sample(rng)
                };
                if let Some(Mul(mul)) = mul {
                    field.extend(quote_spanned! {field.span()=>
                        * #mul
                    });
                }
                if let Some(Div(div)) = div {
                    field.extend(quote_spanned! {field.span()=>
                        / #div
                    });
                }
                if let Some(Add(add)) = add {
                    field.extend(quote_spanned! {field.span()=>
                        + #add
                    });
                }
                if let Some(Sub(sub)) = sub {
                    field.extend(quote_spanned! {field.span()=>
                        - #sub
                    });
                }
                field
            }
            RandomField::DefaultWith(iden) => quote_spanned! {field.span()=>
                #iden
            },
        }
    })
}
