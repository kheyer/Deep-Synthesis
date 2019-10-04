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

WORKDIR /Deep-Synthesis 
RUN pip --no-cache-dir install -r requirements.txt

RUN git clone https://github.com/kheyer/OpenNMT-py 

RUN python build_local/download_model.py

EXPOSE 8501

CMD streamlit run Synthesis/app.py local