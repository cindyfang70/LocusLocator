from flask import Flask, redirect, url_for, request, Response, render_template
from Sequencing2 import *
app = Flask(__name__)

@app.route('/success/<name>')
def success(name):
    return 'welcome %s' % name

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/login', methods=['POST', 'GET'])
def login():
    if request.method == 'POST':
        user = request.form['nm']
        return redirect(url_for('tblastn', name=user))
    else:
        user = request.args.get('nm')
        return redirect(url_for('tblastn', name=user))


@app.route('/tblastn/<name>')
def tblastn(name: str):
    refseqs = blastp(name)
    print("done", flush=True)
    ids = qblast(name)
    get_genomes(ids)
    matches = find_refseq_matches_in_genome(refseqs, "my_genbank.gb")
    results = []
    for match in matches:
        results.append(find_protein_products(find_upstream(convert_to_int(match[1]))))
        results.append(find_protein_products(find_downstream(convert_to_int(match[1]))))
    return render_template("results.html", results=results)




if __name__ == "__main__":
    app.run(debug=True)

