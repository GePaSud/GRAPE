# ==============================
# GRAPE capsule — Docker Full
# ==============================
FROM julia:1.11-bullseye

WORKDIR /grape
COPY . .

RUN apt-get update && apt-get install -y dialog && rm -rf /var/lib/apt/lists/*
RUN julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Ajouter un petit menu shell
COPY start_menu.sh /usr/local/bin/start_menu.sh
RUN chmod +x /usr/local/bin/start_menu.sh

ENTRYPOINT ["/usr/local/bin/start_menu.sh"]
