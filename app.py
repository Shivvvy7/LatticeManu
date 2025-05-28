from flask import Flask, render_template, request
import os
import subprocess
from werkzeug.utils import secure_filename

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
OUTPUT_FOLDER = 'outputs'

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload():
    if 'file' not in request.files:
        return "No file part"
    file = request.files['file']
    if file.filename == '':
        return "No selected file"
    
    filename = secure_filename(file.filename)
    filepath = os.path.join(UPLOAD_FOLDER, filename)
    file.save(filepath)

    # Run MATLAB script and capture output
    try:
        result = subprocess.run([
            "matlab", "-batch",
            f"Main('{filepath.replace('\\', '/')}')"
        ], capture_output=True, text=True)

        print("=== MATLAB STDOUT ===")
        print(result.stdout)
        print("=== MATLAB STDERR ===")
        print(result.stderr)

        # Search for T_avg in MATLAB output
        tavg = None
        for line in result.stdout.splitlines():
            if "Average Temperature of Topmost Layer:" in line:
                tavg = line.split(":")[-1].strip()
                break

        if tavg is None:
            return f"MATLAB ran, but no T_avg found.<br><pre>{result.stdout}</pre><br><pre>{result.stderr}</pre>"

        return render_template("result.html", tavg=tavg)

    except Exception as e:
        return f"Error running MATLAB: {str(e)}"

if __name__ == '__main__':
    app.run(debug=True)
