FROM continuumio/miniconda3

WORKDIR /ICAnnolncRNA

SHELL ["/bin/bash", "--login", "-c"]

COPY environment.yml .
RUN conda env create --file environment.yml

SHELL ["conda", "run", "-n", "lncwebb", "/bin/bash", "-c"]

RUN Rscript -e 'install.packages("LncFinder", repos="https://cloud.r-project.org")'
RUN conda install -y -c conda-forge streamlit==1.26.0
RUN conda install -y -c anaconda pathlib
RUN conda install -y -c anaconda typing
RUN conda install -y -c anaconda nginx
RUN conda install -y -c conda-forge curl

EXPOSE 8501

COPY ICAnnolncRNA .

# Configure nginx reverse proxy for Streamlit app
COPY nginx.conf /etc/nginx/nginx.conf

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "lncwebb", "streamlit", "run", "ICAnnoLncRNA.py"]
