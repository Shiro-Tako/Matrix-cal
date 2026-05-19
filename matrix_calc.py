import numpy as np

EPS = 1e-8


def _to_matrix(values):
    A = np.array(values, dtype=float)
    if A.shape not in {(2, 2), (3, 3)}:
        raise ValueError("Matrix must be 2x2 or 3x3")
    return A


def analyze_matrix(values):
    A = _to_matrix(values)

    eigvals, eigvecs = np.linalg.eig(A)

    real_only = np.all(np.abs(eigvals.imag) < EPS) and np.all(np.abs(eigvecs.imag) < EPS)
    if real_only:
        eigvals = eigvals.real
        eigvecs = eigvecs.real

    rank_v = np.linalg.matrix_rank(eigvecs)
    diagonalizable = rank_v == A.shape[0]

    result = {
        "matrix": A.tolist(),
        "eigenvalues": [complex(v) for v in eigvals],
        "eigenvectors": eigvecs.tolist(),
        "diagonalizable": bool(diagonalizable),
        "reason": None,
        "P": None,
        "P_inv": None,
        "D": None,
    }

    if diagonalizable:
        P = eigvecs
        P_inv = np.linalg.inv(P)
        D = P_inv @ A @ P
        result["P"] = P.tolist()
        result["P_inv"] = P_inv.tolist()
        result["D"] = D.tolist()
    else:
        result["reason"] = "Number of linearly independent eigenvectors is less than matrix size."

    return result


def format_number(x):
    if isinstance(x, complex):
        if abs(x.imag) < EPS:
            return f"{x.real:.6f}"
        sign = "+" if x.imag >= 0 else "-"
        return f"{x.real:.6f} {sign} {abs(x.imag):.6f}i"
    return f"{float(x):.6f}"


def format_matrix(M):
    rows = []
    for row in M:
        rows.append("[ " + ", ".join(format_number(v) for v in row) + " ]")
    return "\n".join(rows)


def build_report(values):
    result = analyze_matrix(values)

    text = "A =\n" + format_matrix(result["matrix"]) + "\n\n"

    text += "1) Eigenvalues\n"
    for i, lmb in enumerate(result["eigenvalues"], start=1):
        text += f"λ{i} = {format_number(lmb)}\n"

    text += "\n2) Eigenvectors (columns of matrix V)\n"
    text += format_matrix(result["eigenvectors"]) + "\n\n"

    text += "3) Diagonalizable?\n"
    text += "Yes\n\n" if result["diagonalizable"] else "No\n\n"

    if result["diagonalizable"]:
        text += "4) P, P^-1, D = P^-1 A P\n"
        text += "P =\n" + format_matrix(result["P"]) + "\n\n"
        text += "P^-1 =\n" + format_matrix(result["P_inv"]) + "\n\n"
        text += "D =\n" + format_matrix(result["D"]) + "\n\n"
        text += "5) If not diagonalizable, reason is shown above.\n"
    else:
        text += "4) Cannot construct P and P^-1 because matrix is not diagonalizable.\n"
        text += f"5) Reason: {result['reason']}\n"

    return text
