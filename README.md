# Sage
Sanger trace alignment

Dependencies
------------

Sage requires Tracy, please intall first:

`https://github.com/gear-genomics/tracy`


Install a local copy for testing
--------------------------------

`git clone https://github.com/gear-genomics/sage.git`

`cd sage`

Setup and run the server
------------------------

The server runs in a terminal

Install the dependencies:

`sudo apt install python python-pip`

`pip install flask flask_cors`

Start the server:

`cd PATH_TO_SAGE/SAGE`

`export PATH=$PATH:/PATH_TO_TRACY/tracy/bin`

`echo $PATH`

`python server/server.py`

Setup and run the client
------------------------

The client requires a different terminal

Install the dependencies:

`cd PATH_TO_SAGE/sage/client`

`sudo apt install npm`

`sudo npm install`

Start the client:

`cd PATH_TO_SAGE/sage/client`

`npm run dev`
