CC = mpicc
CFLAGS = -Wall -Winline -march=native -O3 -fopenmp
DFLAGS = 

INCDIR = inc
OBJDIR = obj
SRCDIR = src

TARGET = GameOfLife
SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJECTS = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SOURCES))

all: clean $(TARGET)

timeit: DFLAGS += -DTIMEIT 
timeit: clean $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) -I$(INCDIR) $(CFLAGS) -o $@ $^ $(DFLAGS)
$(OBJECTS): $(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c -I$(INCDIR) $(CFLAGS) $< -o $@ $(DFLAGS)
clean:
	rm -f $(OBJDIR)/*.o $(TARGET)
