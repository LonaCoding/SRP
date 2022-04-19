import sqlite3 as db


class GeneQuery2:
    def __init__(self):
        self.connection = db.connect('database/GeneQuery2.db')

    def get_clusters(self) -> list:
        cursor = self.connection.cursor()
        find_clusters = cursor.execute('SELECT name from sqlite_master where type="table"')
        clusters_list = [n[0] for n in cursor.fetchall()]
        return clusters_list

    def find_genes(self, cluster: str, gene: str) -> list:
        cursor = self.connection.cursor()
        find = cursor.execute(f'SELECT * from {cluster} WHERE genes = "{gene}"')
        result = cursor.fetchall()
        return result


class FoundGene:
    def __init__(self, row: list):
        self.name = row[0][-1]
        self.p_value = row[0][4]
        self.log2FC = row[0][1]
