# Small-molecule docking image: AutoDock Vina + Meeko (ligand/receptor PDBQT
# prep) + OpenBabel + Biopython, on top of the shared cheminformatics base.
FROM condaforge/miniforge3:latest

WORKDIR /workspace

COPY environment-docking.yml .
RUN mamba env update -n base -f environment-docking.yml && \
    mamba clean -afy
RUN pip install --no-cache-dir meeko

ENV PYTHONUNBUFFERED=1

CMD ["python3"]
