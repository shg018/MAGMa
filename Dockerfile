# LAVA-MS
# Version: 1.0

FROM python:3.8

# Project Enviroment
RUN pip install pipenv
ENV VIRTUAL_ENV=/data/venv
RUN python3 -m virtualenv --python=/usr/bin/python3 $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install Package Libraries
RUN pip install Django django-import-export pandas numpy scipy statsmodels scikit_posthocs openpyxl
# Install R and R Package Libraries
RUN apt-get update && apt-get install -y r-base
RUN R -e "install.packages('limma', repos='http://cran.rstudio.com/')"


# Project Files and Settings
WORKDIR /data/lava
VOLUME ["/data/lava"]
#WORKDIR /workdir/sg2369/LAVA_MS
#VOLUME ["/workdir/sg2369/LAVA_MS"]
#WORKDIR /Users/user/Downloads/TMT_IPMS_Server/src/for_django_3-0
#VOLUME ["/Users/user/Downloads/TMT_IPMS_Server/src/for_django_3-0"]

# Server
EXPOSE 8059
STOPSIGNAL SIGINT
#ENTRYPOINT ["python", "manage.py"]
CMD ["python", "manage.py", "makemigrations"]
CMD ["python", "manage.py", "migrate"]
CMD ["python", "manage.py", "runserver", "0.0.0.0:8059"]