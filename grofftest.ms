.R1
no-label-in-reference
accumulate
.R2
.TL
My first groff doc
.AU
Marius Swane
.AI
NTNU
.AB
This is the abstract. It is a superb summary of all things important!
.AE
.NH
Introduction
.PP
This is a paragraph.
.NH
Next section
.NH 2
This is a subsection
.PP
Here is some more text. I want the word 
.B "bold" " "
in bold text.
Of course you can have
.I italics
as well.
Probably
.UL underscore
as well.
And of course, you can
.BI "combine them!"
.SH
Macros
.IP
This is a weird list
.IP
item
.IP
itme
.IP
item

.de BL
.IP \(bu 2
..
.SH 2
Better looking list
.BL
Now
.BL
the list 
.BL
looks like 
.BL
a list
.PP
A better way to do this is to put all the macros in a macro-file, and source it using .SO.
Preferably at the beginning of the document.
.PP
Let us add some references!
.[
Alesina Artificial states
.]
