FROM docker.io/wcp7un/mmrf_race_docker_r:v4.4.1_d110525

RUN apt-get update

RUN apt-get install -y python3 python3-pip tmux git

#  --break-system-packages not needed until rocker/r-ver:4.4.2
RUN pip install --no-cache-dir radian 

# Other tools
RUN apt-get install -y nano vim neovim htop

RUN apt-get install -y npm
RUN npm install -g tmux-mem tmux-cpu

RUN R -e "install.packages('httpgd')"

# Support for vscode-r in a tmux session
COPY .tmux.conf /root/.tmux.conf

COPY .Rprofile /root/.Rprofile