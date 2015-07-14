# ccmatch [![Build Status](https://travis-ci.org/mkanai/ccmatch.svg?branch=master)](https://travis-ci.org/mkanai/ccmatch)
R package for optimal matching of cases to controls using network flow theory.

## Installation
```{r}
library(devtools)
install_github("mkanai/ccmatch")
```

## Usage
```{r}
> library(ccmatch)
> ccmatch(x.case, x.control, n = 3)
```

## TODO
* Add documentation, example data, detailed usage.
* Using a monotone priority queue [*radix heap*](https://github.com/iwiwi/radix-heap) instead of `std::priority_queue` for perfomance?
  * C++11 is still not supported in every environment yet...

## Reference
* Paul R. Rosenbaum, **Optimal Matching for Observational Studies.** *J. Amer. Statist. Assoc. 84, 408 (Dec., 1989), 1024-1032*

