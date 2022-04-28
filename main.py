from flask import Flask
from flask import render_template
from flask import request
import sqlite3
import os
import database
import template_inserts

app = Flask(__name__)


@app.route('/')
def index():
    with open('/local/www/htdocs/webapp/static/text/main_page/introduction.txt', 'r', encoding='utf-8') as file:
        introduction = file.read()

    return render_template('index.html', introduction=introduction)


@app.route('/pipeline/<int:number>/')
def pipeline(number):
    pipeline_num = number

    # Process relevant text and image files for pipeline page
    template_insert = template_inserts.TemplateInsert(number)
    texts = template_insert.get_text()

    figure_legends = texts[0]
    figure_paragraphs = texts[1]
    methods_results = texts[2]

    figures = template_insert.get_images()
    return render_template('analysis_pipelines.html', pipeline_num=pipeline_num, figure_legends=figure_legends,
                           figure_paragraphs=figure_paragraphs, figures=figures, methods_results=methods_results,
                           zip=zip, type=type, str=str, len=len,
                           list=list)


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

    selected_database = database.GeneQuery(pipeline_num)
    find_gene = selected_database.find_genes(cluster, gene)

    if len(find_gene) == 0:
        found = False
        results = f'{gene} could not be found in {cluster}'
    else:
        found = True
        if pipeline_num == 1:
            results = database.FoundGene1(find_gene)
        else:
            results = database.FoundGene2(find_gene)

    return render_template('genesearch.html', pipeline_num=pipeline_num, gene=gene, cluster=cluster, found=found,
                           results=results, str=str())


@app.route('/references')
def references():
    return render_template('References.html')


if __name__ == '__main__':
    app.run(debug=True)
