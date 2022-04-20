from flask import Flask
from flask import render_template
from flask import request
import sqlite3
import database

app = Flask(__name__)

connection = sqlite3.connect('database/GeneQuery2.db')
cursor = connection.cursor()


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/pipeline/<int:number>/')
def pipeline(number):
    return render_template('analysis_pipelines.html')


@app.route('/pipeline/<int:number>/query', methods=['GET', 'POST'])
def gene_query(number):
    pipeline_num = number
    # TODO: NEED TO RE-STRUCTURE ONCE DATABASE MADE FOR PIPELINE 1
    if pipeline_num == 2:
        database2 = database.GeneQuery2()
        clusters = database2.get_clusters()

    return render_template('genequery.html', pipeline_num=pipeline_num, clusters=clusters)


@app.route('/pipeline/<int:number>/query/search', methods=['POST'])
def search_database(number):
    pipeline_num = number
    cluster = request.form['cluster']
    gene = request.form['gene'].upper()
    # This part might be modified such that gene is a list of genes to search.
    #use loop
    #geneList=list(gene)
    #for g in geneList:
    
    if pipeline_num == 2:
        database2 = database.GeneQuery2()
        find_gene = database2.find_genes(cluster, gene)
        if len(find_gene) == 0: #if search is exhausted yet no match
            found = f'This gene is not present in {cluster}'
        else:
            found = database.FoundGene(find_gene)
            print(found.name)

    return render_template('genesearch.html', pipeline_num=pipeline_num, gene=gene, cluster=cluster, found = found)


@app.route('/references')
def references():
    return render_template('References.html')


if __name__ == '__main__':
    app.run(debug=True)
