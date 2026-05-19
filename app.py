from flask import Flask, render_template_string, request
from matrix_calc import build_report

app = Flask(__name__)

HTML = """
<!doctype html>
<html>
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>Matrix Diagonalization Tool</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 24px; background:#f6f8fc; }
    .card { background:#fff; max-width: 900px; padding: 18px; border-radius: 10px; box-shadow: 0 2px 8px rgba(0,0,0,.08); }
    .grid { display:grid; gap:8px; margin:12px 0; }
    input, select, button { padding:8px; border-radius:8px; border:1px solid #cbd5e1; }
    button { background:#2563eb; color:#fff; border:none; cursor:pointer; }
    pre { background:#0f172a; color:#e2e8f0; padding:14px; border-radius:8px; overflow:auto; }
  </style>
</head>
<body>
  <div class="card">
    <h2>Matrix Assignment Helper</h2>
    <form method="post">
      <label>Size:</label>
      <select name="size" onchange="this.form.submit()">
        <option value="2" {% if size==2 %}selected{% endif %}>2 x 2</option>
        <option value="3" {% if size==3 %}selected{% endif %}>3 x 3</option>
      </select>
      <div class="grid" style="grid-template-columns: repeat({{size}}, 90px);">
        {% for r in range(size) %}
          {% for c in range(size) %}
            <input type="number" step="any" name="a_{{r}}_{{c}}" value="{{matrix[r][c]}}" required>
          {% endfor %}
        {% endfor %}
      </div>
      <button type="submit" name="action" value="compute">Compute</button>
    </form>
    <h3>Result</h3>
    <pre>{{result}}</pre>
  </div>
</body>
</html>
"""


@app.route("/", methods=["GET", "POST"])
def home():
    size = int(request.form.get("size", 2))
    size = 2 if size not in (2, 3) else size

    matrix = []
    for r in range(size):
        row = []
        for c in range(size):
            row.append(request.form.get(f"a_{r}_{c}", "1" if r == c else "0"))
        matrix.append(row)

    result = "Fill matrix values then click Compute."

    if request.method == "POST" and request.form.get("action") == "compute":
        try:
            values = [[float(matrix[r][c]) for c in range(size)] for r in range(size)]
            result = build_report(values)
        except Exception as exc:
            result = f"Error: {exc}"

    return render_template_string(HTML, size=size, matrix=matrix, result=result)


if __name__ == "__main__":
    app.run(debug=True)
