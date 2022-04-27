import os


class TemplateInsert:
    def __init__(self, pipeline_number: int):
        self.text_directory = f'/local/www/htdocs/webapp/static/text/pipeline{pipeline_number}'
        self.images_directory = f'/local/www/htdocs/webapp/static/images/pipeline{pipeline_number}'

    def get_text(self) -> tuple:
        figure_legends_dir = self.text_directory + '/figure_legends'
        figure_paragraphs_dir = self.text_directory + '/figure_paragraphs'
        methods_results_dir = self.text_directory + '/methods_results'

        figure_legends_paths = []
        figure_paragraphs_paths = []
        methods_results = []

        for file in os.listdir(figure_legends_dir):
            path = os.path.join(figure_legends_dir, file)
            figure_legends_paths.append(path)
        for file in os.listdir(figure_paragraphs_dir):
            path = os.path.join(figure_paragraphs_dir, file)
            figure_paragraphs_paths.append(path)
        for file in os.listdir(methods_results_dir):
            path = os.path.join(methods_results_dir, file)
            methods_results.append(path)

        figure_legends = []
        figure_paragraphs = []
        methods_results_text = []

        for text in zip(sorted(figure_legends_paths), sorted(figure_paragraphs_paths)):
            figure_legend = text[0]
            figure_paragraph = text[1]
            with open(figure_legend, 'r', encoding='utf-8') as legend:
                text = legend.read()
                figure_legends.append(text)
            with open(figure_paragraph, 'r', encoding='utf-8') as paragraph:
                text = paragraph.read()
                figure_paragraphs.append(text)

        for text in methods_results_text:
            with open(text, 'r', encoding='utf-8') as file:
                text = file.read()
                methods_results_text.append(text)

        return figure_legends, figure_paragraphs, methods_results_text

    def get_images(self) -> list:
        sorted_figure_paths = sorted(os.listdir(self.images_directory))
        final_paths = []
        for n in range(len(sorted_figure_paths)):
            figure = os.path.join(self.images_directory, sorted_figure_paths[n])
            if figure[-5] == 'a':
                second_figure = os.path.join(self.images_directory, sorted_figure_paths[n + 1])
                final_paths.append([figure, second_figure])
            elif figure[-5] == 'b':
                continue
            else:
                final_paths.append(os.path.join(self.images_directory, sorted_figure_paths[n]))
        return final_paths
