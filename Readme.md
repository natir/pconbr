# Requirements:

- [snakemake](https://snakemake.github.io/)
- [rust toolchains](https://rustup.rs/)
- [jupyter notebook](https://jupyter.org/)
- [conda setup](https://docs.conda.io/en/latest/miniconda.html)

# Reproduce

```
cargo install --git https://github.com/natir/pcon --force
cargo install --git https://github.com/natir/br --force
cargo install --git https://github.com/natir/kmerf --force

cp config.example.yaml config.yaml
# edit config.yaml to set your parameter

jupyter notebook pconbr.ipynb
```

When you are in your jupyter notebook you can run cells in order.
