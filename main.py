from flask import Flask
from flask import render_template
import sqlite3


app = Flask(__name__)

connection = sqlite3.connect('database/GeneQuery2.db')
cursor = connection.cursor()

cursor.execute('SELECT * FROM cluster2 WHERE genes="ETNPPL"')
rows = cursor.fetchall()
print(rows)



@app.route('/')
def index():
    return render_template('index.html')


@app.route('/pipeline/<int:number>/')
def pipeline(number):
    return render_template('analysis_pipelines.html')

@app.route('/pipeline/<int:number>/query')
def gene_query(number):

    return render_template('genequery.html')




@app.route('/references')
def references():
    return render_template('References.html')


if __name__ == '__main__':
    app.run(debug=True)
