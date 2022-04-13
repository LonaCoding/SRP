from flask import Flask
from flask import render_template
import sqlite3
import sqlalchemy
app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')


@app.route('/pipeline/<int:number>/')
def pipeline(number):
    return render_template('analysis_pipelines.html')

@app.route('/pipeline/<int:number>/query')
def gene_query(number):
    db = sqlite3.connect(f'database/GeneQuery{number}.db')
    
    return render_template('index.html')




@app.route('/references')
def references():
    return render_template('References.html')


if __name__ == '__main__':
    app.run(debug=True)
