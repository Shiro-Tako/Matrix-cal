# Matrix-cal

A small web app for the assignment requirements:
- Input a 2x2 or 3x3 matrix
- Compute eigenvalues
- Compute eigenvectors
- Check whether the matrix is diagonalizable
- If diagonalizable, show `P`, `P^-1`, and `D = P^-1 A P`
- If not, show the reason

## Run

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python app.py
```

Open http://127.0.0.1:5000

## Structure

- `app.py`: frontend + web handling
- `matrix_calc.py`: calculation logic (Python)
