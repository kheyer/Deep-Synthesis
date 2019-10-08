FROM continuumio/miniconda3

RUN conda create -n deep_synthesis python=3.6
ENV PATH /opt/conda/envs/deep_synthesis/bin:$PATH
RUN /bin/bash -c "source activate deep_synthesis"
RUN conda install -n deep_synthesis rdkit -c rdkit

RUN git clone https://github.com/kheyer/Deep-Synthesis 

RUN mkdir -p /root/.streamlit

# Streamlit requires something under email in the credentials file
RUN bash -c 'echo -e "\
[general]\n\
email = \"\"\n\
" > /root/.streamlit/credentials.toml'

# CORS currently prevents AWS Load Balancer from
# accessing Streamlit resources
RUN bash -c 'echo -e "\
[global]\n\
logLevel = \"debug\"\n\
\n\
[server]\n\
enableCORS = false\n\
headless = true\n\
" > /root/.streamlit/config.toml'

WORKDIR /Deep-Synthesis 

RUN pip --no-cache-dir install six tqdm==4.30.* \
    future configargparse boto3 pathlib pandas numpy \
    seaborn matplotlib==3.0.0 aiohttp streamlit 

EXPOSE 8501 

CMD streamlit run Synthesis/app.py AWS