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
    if number == 1:
        image = 'statix/image/graoh.png'
        return render_template('analysis_pipelines.html', image)


@app.route('/pipeline/<int:number>/query', methods=['GET', 'POST'])
def gene_query(number):
    pipeline_num = number
    # TODO: NEED TO RE-STRUCTURE ONCE DATABASE MADE FOR PIPELINE 1
    if pipeline_num == 2:
        database2 = database.GeneQuery2()
        clusters = database2.get_clusters()
        gene_count_info = database.GeneCounts()
        genes_per_cluster = gene_count_info.genes_per_cluster
        total_genes = gene_count_info.total_genes

    return render_template('genequery.html', pipeline_num=pipeline_num, clusters=clusters, total_genes=total_genes)


@app.route('/pipeline/<int:number>/query/search', methods=['POST'])
def search_database(number):
    pipeline_num = number
    gene = request.form.get('gene').upper()
    cluster = request.form.get('cluster')
    print(cluster)
    if pipeline_num == 2:
        database2 = database.GeneQuery2()
        find_gene = database2.find_genes(cluster, gene)
        if len(find_gene) == 0:
            found = False
            results = f'{gene} could not be found in {cluster}'


        else:
            found = True
            results = database.FoundGene(find_gene)

    return render_template('genesearch.html', pipeline_num=pipeline_num, gene=gene, cluster=cluster, found=found,
                           results=results)


@app.route('/references')
def references():
    return render_template('References.html')


if __name__ == '__main__':
    app.run(debug=True)
