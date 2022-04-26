import os


class TemplateInsert:
    def __init__(self, pipeline_number: int):
        self.text_directory = f'static/text/pipeline{pipeline_number}'
        self.images_directory = f'static/images/pipeline{pipeline_number}/'

    def get_text(self) -> tuple:
        figure_legends_dir = self.text_directory + '/figure_legends/'
        figure_paragraphs_dir = self.text_directory + '/figure_paragraphs'

        # Obtain list of file paths for text files
        figure_legends_paths = [os.path.join(figure_legends_dir, file) for file in os.listdir(figure_legends_dir)]
        figure_paragraphs_paths = [os.path.join(figure_paragraphs_dir, file) for file in
                                   os.listdir(figure_paragraphs_dir)]

        figure_legends = []
        figure_paragraphs = []
        for text in zip(figure_legends_paths, figure_paragraphs_paths):
            figure_legend = text[0]
            figure_paragraph = text[1]
            with open(figure_legend, 'r') as legend:
                text = legend.read()
                figure_legends.append(text)
            with open(figure_paragraph, 'r') as paragraph:
                text = paragraph.read()
                figure_paragraphs.append(text)
        return figure_legends, figure_paragraphs

    def get_images(self) -> list:
        figures = [os.path.join(self.images_directory, file) for file in os.listdir(self.images_directory)]

        return figures
