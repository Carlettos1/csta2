# Randomizable derive

Fields use existing `Randomizable` implementations unless configured with
`#[csta(range(a..b))]`, `#[csta(len(n))]` for vectors, `#[csta(default)]`,
`#[csta(default = expression)]`, arithmetic modifiers, or `#[csta(after(expression))]`.
`after` first samples that field, then transforms it. Named-field dependencies
are evaluated before consumers; cycles are rejected. Avoid shadowing field names
inside initializer expressions. Dynamic invalid ranges are checked by rand at
sampling time; literal invalid ranges are rejected during compilation.

Enum variants may all carry finite nonnegative numeric `weight` literals, with
at least one positive weight. Zero-weight variants are never sampled. Without
weights, variants are uniform. Unions and empty enums cannot be sampled.

Invalid programs produce compiler diagnostics:

```compile_fail
use csta::csta_derive::Randomizable;
#[derive(Randomizable)]
enum Empty {}
```

```compile_fail
use csta::csta_derive::Randomizable;
#[derive(Randomizable)]
union Unsupported { x: f64 }
```

```compile_fail
use csta::csta_derive::Randomizable;
#[derive(Randomizable)]
enum Zero { #[csta(weight = 0.0)] A }
```

```compile_fail
use csta::csta_derive::Randomizable;
#[derive(Randomizable)]
enum Negative { #[csta(weight = -1.0)] A }
```

```compile_fail
use csta::csta_derive::Randomizable;
#[derive(Randomizable)]
enum Mixed { #[csta(weight = 1.0)] A, B }
```

```compile_fail
use csta::csta_derive::Randomizable;
#[derive(Randomizable)]
struct BadRange { #[csta(range(2..1))] x: f64 }
```

```compile_fail
use csta::csta_derive::Randomizable;
#[derive(Randomizable)]
struct Unknown { #[csta(typo)] x: f64 }
```

```compile_fail
use csta::csta_derive::Randomizable;
#[derive(Randomizable)]
struct Cycle {
    #[csta(default = y)] x: f64,
    #[csta(default = x)] y: f64,
}
```

```compile_fail
use csta::csta_derive::Randomizable;
#[derive(Randomizable)]
struct Missing { #[csta(default = missing_field)] x: f64 }
```
