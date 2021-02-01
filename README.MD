This is a very simple Flask based Progressive Web App that is a demonstration of the different ways to convert a simple Excel spreadsheet calculator tool into a web app 

The live version of this site(*) is deployed using Heroku and available here: https://epidemic-curve-flask-pwa-app.herokuapp.com/

*(Note - if you don't see the Add To Home Screen Button on top left button in your mobile browser - you can add it manually via >Add>Add To Home Screen option)

The actual python file that runs the app is app.py but if you'd like to play around with the calculation functions - have a look at SimpleWebApp.ipynb in this repo. 

The 'original' excel file that it is meant to mimic is located here in OneDrive: https://1drv.ms/x/s!AqZfpTxqClrLhkXfM1yJACw4QZkh?e=DaFtX8 and this folder contains a range of other alternatives explored (e.g. using Excel macros or 3rd party licensing tools to lock/secure excel files vs. this repo which is about converting the spreadsheet into a web app) : https://1drv.ms/u/s!AqZfpTxqClrLhkDA0bcYar6CWQ3b?e=tQ9FYW 

This repo specifically being the example of a DIY build your own web app approach (ref below for other options explored) 

There are 2 other demonstration versions of managed 'web-app' services using other managed web app service tools:
MOLNIFY: https://app.molnify.com/app/epicurveplotter_molnify
GLIDE APPS: https://faint-harbor-0419.glideapp.io/

Read the accompanying Medium Article here: How To Build, Protect & Share A Proprietary Spreadsheet Calculator :
https://zhijingeu.medium.com/how-to-build-protect-share-a-proprietary-spreadsheet-calculator-to-support-lead-generation-da72e7a46cfa

Hit me up on Fiverr if you have a similar challenge and would like some help: https://www.fiverr.com/zhijingeu

Other Notes:

For more on the actual science & math behind the SIR (Susceptible,Infected,Recovered) Model https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model and the Logistic Growth Rate model https://en.wikipedia.org/wiki/Logistic_function#In_medicine:_modeling_of_a_pandemic

Code for the SIR calculation adapted from : Scipython.com https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model

The interactive chart graphics (i.e mouse over tool) were created using the very awesome MPLD3 library in python that translates Matplotlib to D3.JS semi-automatically https://mpld3.github.io/

The code for plotting D3.JS using MPLD3 for this site was adapted from : GitHub MLPD3 Nipunbatra https://github.com/nipunbatra/mpld3-flask

The code that turns this site into a PWA Progressive Web App PWA was based on guidance from PWA Development with Flask https://www.flaskpwa.com/ and https://github.com/umluizlima/flask-pwa




