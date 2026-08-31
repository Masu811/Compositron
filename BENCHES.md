# Compositron - Benchmarks

Disclaimer: The following benchmarks are highly preliminary, as they reflect unoptimized code. Also, keep in mind that Python (NumPy) could be using multi-threading under the hood, so might Julia.

## DBS

| Task                 | Rust    | C++     | Julia   | Python  |
| :---                 | ---:    | :---:   | :---:   | :---:   |
| Ecal-correction      | 155 µs  | 313 µs  | 431 µs  | 1.84 ms |
| BG-subtraction       | 1.53 ms | 1.32 ms | 3.6 ms  | 3.28 ms |
| S-param calculation  | 1.07 µs | 5.89 µs | 4.37 µs | 67.3 µs |

## CDBS

| Task                 | Rust    | C++     | Julia   | Python  |
| :---                 | ---:    | :---:   | :---:   | :---:   |
| Ecal-correction      | 1.96 ms | 3.36 ms | 6.6 ms  | 3.35 ms |
| S-param calculation  | 50.0 µs | 33.0 µs | 106 µs  | 352 µs  |
| Diagonal projection  | 1.84 ms | 1.19 ms | 3.36 ms | 28.3 ms |
| Axis projection      | 1.12 ms | 450 µs  | 90.4 µs | 353 µs  |

## PALS

| Task         | Rust                   | C++                   | Julia               | Python               |
| :---         | ---:                   | :---:                 | :---:               | :---:                |
| 3 LCs, 2 RCs | 102 ms (LM, 140 Iter.) | 266 ms (TR, 24 Iter.) | 73 ms (TR, 68 Iter) | 22 ms (LM, 25 Iter.) |

## Importers

| Format | Rust    | C++     | Julia   | Python  |
| :---   | ---:    | :---:   | :---:   | :---:   |
| SLOPE  | 5.06 ms | 4.73 ms | 16.5 ms | 11.8 ms |
