#!/usr/bin/env python3

'''
NAME:	    restriction_GUI.py
Purpose:	This module can be used to display an restriction analyses
Author:	 	Ronald Wedema
Created:	22/03/2012
Copyright:  (c) Ronald Wedema 2012
License:	GNU General Pulic License (GPL) version 3
version:	2.0
'''

from tkinter import *
from tkinter.filedialog import askopenfilename
from math import log

#import hier je API performRestriction
#from t3_import performRestriction, jouw script moet dus t3.py heten en er moet een functie zijn die performRestriction heet die al het werk doet
#zie functie digest over wat de argumenten zijn van performRestriction en wat de return is...
from t3 import performRestriction


#RadioButtons
MODES = { "EcorI" : ["GAATTC", 1], "BamHI" : ["GGATCC", 1], "EcoRII" : ["CCWGG", 0], "HindIII":  ["AAGCTT", 1], "TaqI" : ["TCGA", 1], "NotI" : ["GCGGCCGC", 2], "HinfI" : ["GANTCA", 1], "Sau3A" : ["GATC", 0], "HaeIII" : ["GGCC", 2], "EcoRV" : ["GATATC", 3], "PstI" : ["CTGCAG", 4], "XbaI" : ["TCTAGA", 1], "MboI" : ["GATC", 0], "mst2" : ["CCTNAGG",2]}

#ambigious nucleotides!! dus regex gebruiken, zie bv. EcoRII en HinfI
#N = C or G or T or A
#W = A or T


enzymes = sorted(MODES.keys())
test_fragments = {"frag1": [2,3,4], "frag2" : [3,4,5], "frag3" : [4,5,9]}

button_num = 0
button_dic = {}
var_list = []
cb_list = []
DNAFile = ""

def main():
	root = Tk()

	#RootFrame
	frame = Frame(root)
	frame.pack()

	buildMenus(frame)

	root.mainloop()


def buildMenus(frame):
	#Topmenu
	menubar = Frame(frame,relief=RAISED,borderwidth=1)
	menubar.pack(side=TOP)

	#RadioButtonsFrame
	buttonFrame = Frame(frame,relief=RAISED,borderwidth=1)
	buttonFrame.pack(side=LEFT)

	#TextFrame
	textFrame = Text(frame,relief=RAISED,borderwidth=1,width=80, height=100)
	textFrame.pack(side=LEFT)

	#CanvasFrame
	#Adapt this size for the canvas
	canvasFrame = Canvas(frame,background = "black",relief=RAISED,borderwidth=1, width=800, height=1000)
	canvasFrame.pack(side=LEFT)

	#DigestButtonFrame
	ExecuteFrame = Frame(frame,relief=RAISED,borderwidth=1)
	ExecuteFrame.pack(side=BOTTOM)

	#TopMenu Buttons
	mb_file = Menubutton(menubar,text='file')
	mb_file.pack(side=LEFT)
	mb_file.menu = Menu(mb_file)

	mb_file.menu.add_command(label='open', command=lambda: openFile(textFrame))
	mb_file.menu.add_command(label='close', command=frame.quit)

	mb_edit = Menubutton(menubar,text='edit')
	mb_edit.pack(side=LEFT)
	mb_edit.menu = Menu(mb_edit)
	mb_edit.menu.add_command(label='copy')
	mb_edit.menu.add_command(label='paste')

	mb_help = Menubutton(menubar,text='help')
	mb_help.pack(padx=25,side=RIGHT)

	mb_file['menu'] = mb_file.menu
	mb_edit['menu'] = mb_edit.menu

	#DigestButton
	ExecuteButton = Button(ExecuteFrame, text="Digest!", command=lambda: digest(canvasFrame))
	ExecuteButton.pack(side=RIGHT)

	#ClearGelButton
	ClearGelButton = Button(ExecuteFrame, text="Clear Gel!", command=lambda: reset(canvasFrame))
	ClearGelButton.pack(side=RIGHT)
	
	#ClearTextButton
	ClearTextButton = Button(ExecuteFrame, text="Clear text!", command=lambda: clearText(textFrame))
	ClearTextButton.pack(side=RIGHT)

	buttons(frame)

def reset(frame):
	frame.delete(ALL) #remove all items
	
	
def clearText(frame):
	frame.delete(1.0, END) #remove all Text


def openFile(textFrame):
	'''gebruikt de GUI om een bestand naam te vragen en retourneert deze'''

	global DNAFile
	filename = askopenfilename(filetypes=[("allfiles","*")])
	textFrame.insert(INSERT, filename + "\n")

	readFasta(filename, textFrame)
	DNAFile = filename


def readFasta(filename, textFrame):
	'''gebruikt de geselecteerde file om de inhoud te printen in de GUI'''

	for line in open(filename, 'r').readlines():
		textFrame.insert(INSERT, line)
		textFrame.pack()


def digest(canvasFrame):
	'''deze functie roept de API aan en wordt direct gestart na het klikken op de digest button'''

	#haal de enzymen op die geselecteerd zijn
	enzymes_selected = finished()

	#roep de API aan, geef path DNA file en lijst met geselecteerde enzymen mee
	result_dict = performRestriction(DNAFile, enzymes_selected)

	#teken de fragmenten op het canvas
	#disable this function call when no test fragments should be used!!
	drawgel(canvasFrame, test_fragments)

	#enable below function when your own fragments should be drawn!!
	drawgel(canvasFrame, result_dict)


def calcDistance(mw):
	distance = 1 / (math.log(mw))
	return distance


def drawgel(canvasFrame, fragments):
	'''deze functie tekend de verschillende monsters en hun fragmenten op het canvas'''

	#Line: x0, y0, x1, y1
	#frame = x200, y400 (total dimension canvas, adapt to your needs)

	#determine total gel lanes ( =monsters )
	lanes = len(fragments)
	xshift = 800 / lanes

	x_text = 50
	y_text = 10

	Frag_max = 0

	#get max fragment size
	for key, value in fragments.items():		
		if (max(value) > Frag_max):
			Frag_max = max(value)


	yshift = 900 / Frag_max

	x0 = 0
	x1 = 100
	y0 = 20
	y1 = y0

	coord = x0, y0 , x1, y1

	#place monster names above the lanes
	for key, value in fragments.items():
		
		text_pos = x_text, y_text
		canvasFrame.create_text(text_pos, fill= "yellow", text=key)

		x_text = x_text + xshift

		#draw each fragment
		for fragment in value:
			y0 = 20
			y1 = y0
			y0 = y0 + (fragment * yshift)
			y1 = y0

			coord = x0, y0 , x1, y1

			canvasFrame.create_line(coord, fill='red', width=3)
			canvasFrame.pack()

		x0 = x0 + xshift
		x1 = x1 + xshift

def buttons(frame):
	""" one checkbutton for each item in the list """

	for enzyme in enzymes:
		button_dic[button_num] = enzyme
		var = IntVar()
		var.set(0) ## off=default
		var_list.append(var)

		#create button
		b = Checkbutton(frame, text = enzyme, variable=var)
		b.pack()

		#add to list
		cb_list.append(b)


def finished():
	""" called from the print button"""
	enzymes_selected = []

	ctr = 0

	#check if button has been selected
	for var in var_list:
		v = var.get()

		if v:
			#what to return only enzyme name
			enzymes_selected.append(enzymes[ctr])

		ctr += 1

	return enzymes_selected

if __name__ == "__main__":
	main()
