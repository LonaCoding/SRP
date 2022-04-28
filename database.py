import sqlite3 as db

''' Class to interact with both database files'''


class GeneQuery:
    def __init__(self, database_number: int):  # Depending on the url the user is on, the class is initialised with 1
        # or 2 representing the pipeline.
        self.database_number = database_number
        self.connection = db.connect(f'/local/www/htdocs/webapp/database/GeneQuery{self.database_number}.db')  # This
        # number is used to connect with a specific database

    # A method that obtains the clusters from the database
    def get_clusters(self) -> list:
        cursor = self.connection.cursor()
        find_clusters = cursor.execute('SELECT name from sqlite_master where type="table"')
        clusters_list = [n[0] for n in cursor.fetchall()]
        return clusters_list

    # A method that tries to find a selected gene from a selected cluster (cluster and gene are selected from the
    # genequery webpage through a POST request)
    def find_genes(self, cluster: str, gene: str) -> list:
        cursor = self.connection.cursor()

        find = cursor.execute(f'SELECT * from {cluster} WHERE gene = "{gene}"')
        result = cursor.fetchall()
        return result

    # Method that counts all the genes within a cluster
    def gene_counts(self, cluster: str):
        cursor = self.connection.cursor()
        find = cursor.execute(f'SELECT COUNT(*) FROM {cluster}')
        result = cursor.fetchall()
        return result


''' A class to set attributes for a given gene and cluster along with its experimental data for database 2'''


class FoundGene2:
    def __init__(self, row: list):
        self.name = row[0][0]
        self.p_value = row[0][-1]
        self.log2FC = row[0][2]


''' A class to set attributes for a given gene and cluster along with its experimental data for database 1'''


class FoundGene1:
    def __init__(self, row: list):
        self.name = row[0][0]
        self.p_value = row[0][-1]
        self.log2FC = row[0][2]


''' A class that couts all the genes detected from all cells '''


class GeneCounts(GeneQuery):
    def __init__(self, database_number):
        super().__init__(database_number)
        self.genes_per_cluster = [self.gene_counts(n)[0][0] for n in self.get_clusters()]
        self.total_genes = sum(self.genes_per_cluster)
