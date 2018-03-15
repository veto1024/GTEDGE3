#!/usr/bin/python

############################################3
#
#
#    Graphical interface files
#
#
############################################

import Tkinter as tk

NORM_FONT=("Helvetica",10)

def trashCont():
    return
def popup(msg,title=""):
    popup=tk.Tk()
    popup.wm_title(title)
    label=tk.Label(popup,text=msg,font=NORM_FONT)
    label.pack(side="top",fill="x",pady=10)
    B1=tk.Button(popup,text="Okay",command=popup.destroy)
    B1.pack()
    popup.mainloop()

if __name__=="__main__":
    
    pass