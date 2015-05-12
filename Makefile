#
#  Makefile for MAC OS (cocoa)
#
#     make all   (construction de l'executable)
#     make clean (effacement des fichiers objets et de l'executable)
#
#  A adapter en fonction des ordinateurs/environnements 
#  Compilateur, edition de liens, 
#
GLFWINCLUDE  = /Users/Ed/Documents/EPL/C/Devoir3/glfw-3.1/include
GLFWLIBRARY  = /Users/Ed/Documents/EPL/C/Devoir3/glfw-3.1/src/libglfw3.a
#
CC       = gcc  
LD       = gcc
CFLAGS   = -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk -mmacosx-version-min=10.6 -I$(GLFWINCLUDE) -O3 -Dgraphic -Wall -g
LFLAGS   = -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk -mmacosx-version-min=10.6 -Wl,-search_paths_first -Wl,-headerpad_max_install_names -g -Wall -O3 -framework AGL -framework Cocoa -framework OpenGL -framework IOKit -framework CoreFoundation -framework CoreVideo
LIBS     = $(GLFWLIBRARY) -lm 
#
#
PROG     = myFem
LISTEOBJ = \
    main.o
# ATTENTION... aucun caractere apres le caractere de continuation "\"
#
# compilation
#
.c.o :
	$(CC) -c  $(CFLAGS) -o $@ $<
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# dependances
#
all        : $(PROG)
cube.o     : cube.c 
main.o     : main.c cube.c
#
# edition de lien
#
$(PROG) : $(LISTEOBJ)
	$(LD) -o $(PROG) $(LFLAGS) $(LISTEOBJ) $(LIBS)
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# effacement des fichiers intermediaires
#
clean :
	rm -vf $(PROG) $(LISTEOBJ) core a.out
# ATTENTION... la ligne precedente doit commencer par une tabulation
#
# ATTENTION... il faut une ligne vide a la fin du fichier.


