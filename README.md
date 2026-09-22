# symcas

Fast, deterministic computer algebra system in Rust.

`symcas` provides a canonical expression arena with hash-consing, exact arithmetic over arbitrary precision rational numbers, ordered sparse multivariate polynomials, symbolic differentiation, series expansions, directed simplifications, roundtrip plain parsing, and LaTeX rendering.

## Key Features

- **Canonical Construction & O(1) Equality**: Expressions are automatically normalized upon construction through flattening, term sorting, exact arithmetic folding, and like-term combination. Equivalent expressions within the same `Context` share an identical arena handle (`Expr::raw_id`).
- **Exact Numeric Domain**: Arbitrary-precision integers and rationals backed by `num-bigint` / `num-integer` with inline `i64` optimization for small integers.
- **Multivariate Polynomials**: Ordered sparse distributed polynomials generic over coefficient rings, supporting addition, multiplication, exact division, polynomial remainder sequences (PRS), multivariate GCD, and square-free factorization.
- **Symbolic Calculus**: Symbolic differentiation over elementary functions (`sin`, `cos`, `tan`, `exp`, `log`, `sqrt`, `abs`), Taylor series expansions, algebraic expansion with size bounds, and cancellation.
- **Assumptions & Predicates**: Three-valued logic (`True`, `False`, `Unknown`) with deduplicated predicates (`Positive`, `NonZero`, `Integer`, etc.).
- **Deterministic Roundtrip**: Machine-readable `plain` formatting that parses back to the exact same expression node (`parse(plain(e)) == e`), along with publication-ready `latex` output.

## Architecture

The system is organized into modular crates:

| Crate | Role |
|---|---|
| `symcas` | High-level facade crate: prelude, `parse`, `plain`, and `latex` APIs |
| `cas-domain` | Numeric core: `Integer`, `Rational`, `Ring`, and `Field` traits |
| `cas-poly` | Multivariate polynomial algebra, ordering (`Lex`, `DegRevLex`), and GCD |
| `cas-expr` | Expression arena, DAG canonicalization, calculus, and assumptions |
| `cas-print` | Deterministic `plain` and `latex` printers |
| `cas-parse` | Recursive descent parser for `plain` mathematical syntax |

## Quickstart

Add `symcas` to your `Cargo.toml`:

```toml
[dependencies]
symcas = "0.1.0"
```

### 1. Expressions & Canonical Form

```rust
use symcas::prelude::*;

let ctx = Context::new();
let x = ctx.sym("x");
let y = ctx.sym("y");

// Arithmetic operations construct canonical forms on the fly
let e = (x.clone() + y.clone()).pow(3) - x.clone().pow(3) - y.clone().pow(3);
assert_eq!(symcas::plain(&ctx, &e), "(x + y)^3 - x^3 - y^3");

// Hash-consing: expressions with identical canonical structure share identical node IDs
let a = x.clone() + x.clone();
let b = ctx.int(2) * x.clone();
assert!(a == b);
```

### 2. Parsing & Roundtrip Guarantee

```rust
use symcas::prelude::*;

let ctx = Context::new();
let e = symcas::parse(&ctx, "x + x + 1/2*x").unwrap();
assert_eq!(symcas::plain(&ctx, &e), "5/2*x");

// Parsing plain output yields the exact same canonical node
let back = symcas::parse(&ctx, "5/2*x").unwrap();
assert!(back == e);
```

### 3. Expansion, Simplification & Substitution

```rust
use symcas::prelude::*;

let ctx = Context::new();
let (x, y) = sym!(&ctx, x, y);

// Expand products and powers
let e = (x.clone() + y.clone()).pow(2);
let expanded = ctx.expand(&e);
assert_eq!(symcas::plain(&ctx, &expanded), "x^2 + y^2 + 2*x*y");

// Algebraic substitution
let sub = ctx.subst(&expanded, &[("x", ctx.int(1))]);
assert_eq!(symcas::plain(&ctx, &sub), "1 + y^2 + 2*y");

// Trigonometric simplification
let trig = ctx.call("sin", &[x.clone()]).pow(2) + ctx.call("cos", &[x.clone()]).pow(2);
let sim = ctx.simplify(&trig);
assert_eq!(symcas::plain(&ctx, &sim), "1");
```

### 4. Differentiation & Taylor Expansion

```rust
use symcas::prelude::*;

let ctx = Context::new();
let x = ctx.sym("x");

// Symbolic derivative d/dx (sin(x) * x^2)
let f = ctx.call("sin", &[x.clone()]) * x.clone().pow(2);
let df = ctx.diff(&f, &x);
assert_eq!(symcas::plain(&ctx, &df), "2*x*sin(x) + x^2*cos(x)");

// Taylor series around x = 0 to order 4 for exp(x)
let ex = ctx.call("exp", &[x.clone()]);
let s = ctx.taylor(&ex, &x, 0, 4);
assert_eq!(symcas::plain(&ctx, &s), "1 + x + 1/24*x^4 + 1/6*x^3 + 1/2*x^2");
```

### 5. Rational Reduction & Polynomial Operations

```rust
use symcas::prelude::*;

let ctx = Context::new();
let x = ctx.sym("x");

// Rational cancellation via polynomial GCD
let frac = (x.clone().pow(2) - ctx.int(1)) / (x.clone() - ctx.int(1));
let reduced = ctx.cancel(&frac);
assert_eq!(symcas::plain(&ctx, &reduced), "1 + x");

// Factoring polynomials
let poly = x.clone().pow(2) - ctx.int(4);
let factored = ctx.factor(&poly);
assert_eq!(symcas::plain(&ctx, &factored), "(-2 + x)*(2 + x)");
```

## License

Licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT License ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.
