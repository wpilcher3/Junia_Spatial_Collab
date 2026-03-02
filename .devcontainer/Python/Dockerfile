FROM docker.io/wcp7un/mmrf_race_docker_python:v3.13_d120225b

RUN apt-get update
RUN apt-get install -y tmux git nano vim neovim htop npm && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Other tools
RUN npm install -g tmux-mem tmux-cpu


# Support for vscode-r in a tmux session
COPY .tmux.conf /root/.tmux.conf
