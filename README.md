# cl-ceigen-lite

[C2FFI](https://github.com/rpav/c2ffi/) / [cl-autowrap](https://github.com/rpav/cl-autowrap) based wrapper for [ceigen_lite](https://github.com/digikar99/ceigen_lite).

This project has been to put to use in [numericals](https://github.com/digikar99/numericals).

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [ceigen](#ceigen)
- [ceigen-rand](#ceigen-rand)
    - [Set seed](#set-seed)
    - [Floating point distributions](#floating-point-distributions)
    - [Integer distributions](#integer-distributions)

<!-- markdown-toc end -->


## ceigen

Undocumented. See [ceigen\_lite/ceigen\_lite.h](./ceigen_lite/ceigen_lite.h) for a list of functions.

## ceigen-rand

Lightweight wrapper around random number generators of ceigen-lite. These functions are not thread safe. Use reentrant `*-r` versions of functions in ceigen-lite to write your own. Refer to [EigenRand](https://bab2min.github.io/eigenrand/v0.6.0.rc2/en/list_of_supported_distribution.html) for documentation.

### Set seed

- seed

### Floating point distributions

- beta
- beta!
- cauchy
- cauchy!
- chi-squared
- chi-squared!
- exponential
- exponential!
- extreme-value
- extreme-value!
- fisher-f
- fisher-f!
- gamma
- gamma!
- lognormal
- lognormal!
- normal
- normal!

### Integer distributions

These only support (signed-byte 32) given EigenRand's limitations.

- bernoulli
- bernoulli!
- binomial
- binomial!
- geometric
- geometric!
- negative-binomial
- negative-binomial!
- poisson
- poisson!
- uniform-int
- uniform-int!
