from flask import Flask
from flask import render_template
from flask import request
import sqlite3
import database

app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/pipeline/<int:number>/')
def pipeline(number):
    return render_template('analysis_pipelines.html')


@app.route('/pipeline/<int:number>/query', methods=['GET', 'POST'])
def gene_query(number):
    pipeline_num = number

    selected_database = database.GeneQuery(pipeline_num)
    clusters = sorted(selected_database.get_clusters())
    gene_count_info = database.GeneCounts(pipeline_num)
    total_genes = gene_count_info.total_genes

    return render_template('genequery.html', pipeline_num=pipeline_num, clusters=clusters, total_genes=total_genes,
                           len=len)


@app.route('/pipeline/<int:number>/query/search', methods=['POST'])
def search_database(number):
    pipeline_num = number
    gene = request.form.get('gene').upper()
    cluster = request.form.get('cluster')
    print(cluster)
    if pipeline_num == 2:
        database2 = database.GeneQuery(pipeline_num)
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
