# GOMonteCarlo
This software package simulates the reduction of graphene oxide and calculates the optical response of the generated sp2 domains

Prerequisites:
GTK+3
gnuplot


Installation Instructions:

Linux:
compile with this command
gcc GOMonteCarlo.c -o GOMonteCarlo `pkg-config --cflags --libs gtk+-3.0` -lm -lpthread 

To Run:
./GOMonteCarlo


Mac:
Use this script to install the prerequisites:
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install gtk+3
install gnuplot
BASEDIR="$( dirname "$0" )"
cd "$BASEDIR"

Then compile
gcc GOMonteCarlo.c -o GOMonteCarlo `pkg-config --cflags --libs gtk+-3.0` -lm -lpthread 
./GOMonteCarlo

Windows:
At present, the program is missing some precompiler flags to switch threading to windows method of threading.
Also, there are some issues with using GTK+ in windows, but I have compiled GTK+ code in Windows before.
If you would like to run this program on a Windows machine please contact me.
