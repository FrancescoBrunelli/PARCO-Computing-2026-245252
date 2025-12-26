# Results
This directory contains all the performance metrics collected from the experiments, including **execution time**, **perf**, and **cachegrind** data.  
The results are organized by **matrix**, **scheduling policy**, and **memory configuration**.

---

## Directory Structure

``` text
results
├── matrix_name         # One for each matrix
│   ├── dynamic
│   │   ├── 128GB
│   │   │   ├── perf
│   │   │   └── runtime
│   │   ├── 32GB
│   │   │   ├── perf
│   │   │   └── runtime
│   │   └── 64GB
│   │       ├── perf
│   │       └── runtime
│   ├── guided
│   │   ├── 128GB
│   │   │   ├── perf
│   │   │   └── runtime
│   │   ├── 32GB
│   │   │   ├── perf
│   │   │   └── runtime
│   │   └── 64GB
│   │       ├── perf
│   │       └── runtime
│   └── static
│       ├── 128GB
│       │   ├── perf
│       │   └── runtime
│       ├── 32GB
│       │   ├── perf
│       │   └── runtime
│       └── 64GB
│           ├── perf
│           └── runtime
└── cachegrind_logs
    └── matrix_name         # One for each matrix
        ├── dynamic
        ├── guided
        └── static
```