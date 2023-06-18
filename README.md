# solpos

Solar Position Algorithm ([NREL SPA](https://midcdmz.nrel.gov/spa/)) implemented in Rust, with Python bindings.

## Setup & Installation

### Prerequisites

1. [Rust](https://www.rust-lang.org/tools/install)
2. [Python 3.10](https://www.python.org/downloads/release/python-3100/)
3. [Poetry](https://python-poetry.org/docs/#installation)
4. [Maturin](https://github.com/PyO3/maturin#readme)
5. A Python virtual environment

### Steps

1. Clone the repository:

```
git clone https://github.com/username/solpos.git
cd solpos
```

2. Setup your Python virtual environment:

```
python3 -m venv venv
source venv/bin/activate
```

3. Install the project dependencies using Poetry:

```
poetry install

```

4. To build the project, run:

```
maturin develop
```

5. To run the Rust tests:

```
cargo test
```

6. To run the Python tests:

```
pytest
```

7. To run the Rust benchmarks:

```
cargo bench

```
