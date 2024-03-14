# Antrag

This repository contains an implementation of the [Antrag][antrag-eprint]
lattice-based signature scheme (as presented at [ASIACRYPT 2023][ac2023])
with several additional optimization discussed in particular in the
presentation at the [NIST 5th PQC Standardization
Conference][nist-5thpqc].

The implementation is based on the [earlier implementation of
Antrag][antrag-implorig] and on the Round 3 implementation of
[Falcon][falcon-impl], with extensive modifications to support more
flexible parameter selection. It is mostly the work of Jade Guiton (with
contributions by Mehdi Tibouchi), carried out as part of her internship
at the NTT Social Informatics Laboratories.

[antrag-eprint]: https://eprint.iacr.org/2023/1335
[ac2023]: https://asiacrypt.iacr.org/2023/
[nist-5thpqc]: https://csrc.nist.gov/events/2024/fifth-pqc-standardization-conference
[antrag-implorig]: https://github.com/mti/antrag
[falcon-impl]: https://falcon-sign.info/impl/falcon.h.html

## Directory structure

- `antrag`: Code from the original Antrag implementation (modified)
- `falcon`: Code from Falcon (modified)
- `ntru`: Code for new NTRUSolve implementation
- `supercop`: Code for integration with SUPERCOP
- `tests`: Code for tests and benchmarks
- `scripts`: SageMath scripts for testing and generating lookup tables

The following 3 directories are automatically created as part of the build
process.
- `gen`: Headers generated from schema parameters by `scripts/gen_headers.sage`
- `build`: Intermediary compilation artifacts
- `supercop-build`: Generated code for inclusion in SUPERCOP

## Makefile targets

- `make`: Builds the main end-to-end demo program / benchmark
- `make antrag1f`: Same as above, selecting the `antrag1f` parameter set
  (the full list of parameter sets that can be chosen is as follows:
  `antrag1f`, `antrag1s`, `antrag2f`, `antrag2s`, `antrag3f`, `antrag3s`,
  `antrag4f`, `antrag4s`, `antrag5f`, `antrag5s`)
- `make test-main`: Like plain `make`, but also runs the `main` program
- `make test-ntt`: Runs a unit test of the NTT routines
- `make test-poly`: Runs a unit test of the FFT and complex polynomial manipulation routines
- `make test-codec`: Runs an end-to-end test including key and signature encoding/decoding
- `make test-rns`: Runs a unit test of the RNS and NTT-under-RNS routines
- `make benchmark`: Runs a benchmark for various Antrag routines
- `make build-supercop`: Flattens the source tree into `supercop-build` for use in SUPERCOP (see below)
- `make clean`: Reset build tree

## How to benchmark with SUPERCOP

- Download the SUPERCOP archive, extract it in a new folder
- Remove all lines but the first in `okcompilers/c` (to avoid trying many compiler configurations)
- Create the directory `crypto_sign/antrag`
- Run `./do-part init`
- Run `./do-part gmp`
- Run `./do-part crypto_stream chacha20`
- Run `./do-part crypto_rng chacha20`

The final three steps must be redone if you change the code in Antrag.
- Run `make build-supercop` in the Antrag tree
- Symlink or copy the generated directory (`supercop-build`) into the SUPERCOP tree at `crypto_sign/antrag/ref`
(If you choose symlinking, this step only has to be done once.)
- Run `./do-part cypto_sign antrag`

The resulting benchmark database file should be at `bench/<folder>/data`.
