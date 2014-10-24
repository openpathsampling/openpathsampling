
 /*                 xising
 
   Xwindow program to simulate the Ising model using microcanonical 
   algorithms.  Many "demons" circulate around the lattice
   trying to flip spins. 

  Michael Creutz   
  creutz@wind.phy.bnl.gov 

  this is version 2.2 April 29, 1994
  version 1.0 was for SPARC's only
  version 2.0 for generic X systems
  version 2.1 colored the buttons and thermometer fluid
  version 2.2 fixes bit ordering problems for DEC alpha's

compiling:  cc -O -o xising xising.c -lm -lX11

If you find some machine supporting X that it does not work on, 
please let me know.

I have not included a canonical algorithm because 1) for the
local algorithm it is noticeably slower, and 2) I don't believe
that one would visually be able to distinguish any other difference 
beyond a non-fluctuating thermometer. */

# include <X11/Xlib.h>
# include <X11/Xutil.h>
# include <X11/Xos.h>
# include <X11/Xatom.h> 
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <limits.h>

 /* lattice dimensions:
        NCOLS must be a multiple of WORDSIZE 
        NROWS should not be a multiple of 29 */
# define NROWS 256
# define NCOLS 256

# define WORDSIZE  (8*sizeof(long)) 
# define HWORDS (NCOLS/WORDSIZE)
# define VOLUME (NROWS*NCOLS)
# define NDEMONS (2*NROWS*HWORDS)
# define LEFTBIT  (((long) 1)<<(WORDSIZE-1))

long field[NROWS][HWORDS];      /* to store the system */

 /* some convenient macros */
# define nextrow (row+1-NROWS*(row>=NROWS-1))
# define nextcol (col+nshift*(1-HWORDS*(col==topcol)))
# define previousrow (row-1+NROWS*(row<1))
# define previouscol (col-nshift*(1-HWORDS*(col==bottomcol)))

long mrand48(),lrand48(); /* generates a random word, lrand is non-negative */
double drand48(); /* generates a random double */
char stringbuffer[100];
static char *progname;
char bits[256]; /* for counting bits */

/* bitcount counts the set bits in a word; call setup() before using */
# define bitcount8(i)  (bits[(i)&255])
# define bitcount16(i) (bitcount8(i)+bitcount8((i)>>8))
/* check if long int is 32 bits */
#if ((0x7FFFFFFF)==(LONG_MAX))  /* 32 bits */
# define bitcount(i)   (bitcount16(i)+bitcount16((i)>>16))
# define RND mrand48()
#else  /* assume 64 bits */
# define bitcount32(i)   (bitcount16(i)+bitcount16((i)>>16))
# define bitcount(i)   (bitcount32(i)+bitcount32((i)>>32))
# define RND (mrand48()^(mrand48()<<32))
#endif

/* things to monitor temperature, etc. */
double beta; /* the inverse temperature */
long energycount=0,volumecount=0;
int heater=0,cooler=0,aboutmapped=0;
int boundary=0,algorithm=0, paused=0, iteration=0;
  /* variables to allow for peculiar bit orderings on some machines */
int nshift,topcol,bottomcol;

/* for local algorithm; use two demon bits */
long demon0=0,demon1=0,checkermask;
long work0[HWORDS],work1[HWORDS];

/* stuff for the cluster algorithm; use three demon bits */
# define DEMONBITS 3
long cluster[NROWS][HWORDS];    /* to store cluster shape */
long sadx[NROWS][HWORDS],sady[NROWS][HWORDS]; /* labels sad demons */
int unchecked[NROWS][HWORDS]; /* to reduce redundant calculations */
int boxtop,boxbottom,boxleft,boxright; /* will enclose cluster */
long demon[DEMONBITS][NDEMONS]; /* for demons in cluster algorithm */
long newdemon[DEMONBITS][NDEMONS];
long demonindex=0; /* for shuffling demons */
long clustergrowth=0;
int direction=1; /* for vertical sweeps through cluster */

/* various window stuff */
Display *display;
int screen;
Window window,quitbutton,pausebutton,algobutton,heatwindow,boundwindow,about;
GC gc, gcr, gcxor, gccolor;
XImage *spinimage,*clusterimage;
XFontStruct *font=NULL;
int font_height;
unsigned int windowwidth, windowheight;
XSizeHints size_hints;

/* dimensions for control boxes */
#define HEATERHEIGHT 60
#define HEATERWIDTH 96
#define BOUNDHEIGHT 80
#define BOUNDWIDTH  110

main(argc,argv)
int argc;
char **argv;
{unsigned int width, height;
 int row,col;
 XEvent report;
 progname=argv[0];
 openwindow(argc,argv);
/*
 setup();
  
 for (col=0;col<HWORDS;col++)
 for (row=0;row<NROWS;row++)
  {field[row][col]=~(RND*(drand48()>.77));
   cluster[row][col]=0;
   unchecked[row][col]=1;
  }
 beta=0.5*log(1.0+sqrt(2.0)); 
 refreshdemons();
*/
/* loop forever, looking for events */
 while(1)
   {/*if (clustergrowth) 
      while (clustergrowth) growcluster(); 
   else */
     XClearArea(display, window, 0,0, windowwidth, windowheight,0);
     startcluster();
     
         
   if (paused|XPending(display))
    {XNextEvent(display,&report); 
     switch (report.type)
      {case Expose:
        if (report.xexpose.count!=0) break; /* more in queue, wait for them */
/*        repaint();  */
        break;
       case ConfigureNotify:
        width=report.xconfigure.width;
        height=report.xconfigure.height;
        if ((width<size_hints.min_width)||(height<size_hints.min_height))
            {fprintf(stderr,"%s: window too small to proceed.\n",progname);
             XUnloadFont(display,font->fid);
             XFreeGC(display,gc); 
             XCloseDisplay(display);
             exit(1);
            } 
        break; 
/*       case ButtonPress:
        if (report.xbutton.window==quitbutton)
            {XUnloadFont(display,font->fid);
             XFreeGC(display,gc); 
             XCloseDisplay(display);
             exit(1);
            } 
        else if (report.xbutton.window==pausebutton)
            {paused=1-paused;
             mypause();
            }
        else if (report.xbutton.window==algobutton)
            {algorithm=1-algorithm;
             setalg();
             myalgo();
            }
        else if (report.xbutton.window==heatwindow)
            {setheater((report.xbutton.y*3)/HEATERHEIGHT);
            }
        else if (report.xbutton.window==boundwindow)
            {setboundary((report.xbutton.y*4)/BOUNDHEIGHT);
            }
        else if (0==aboutmapped)
            {XMapWindow(display,about);
             aboutmapped=1;
            }
        else if (1==aboutmapped)
            {XUnmapWindow(display,about);
             aboutmapped=0;
            }
        break;*/
       default:
        break;
      }
    }
  }
}

startcluster()
{
  int x1,x2,y1,y2,i;
  
  y1=x1=0;
  for(i=0;i<30;i++) {
    x2 =random()*256;
    y2 =random()*256;
    XDrawLine(display,window,gc,x1,y1,x2,y2);
    x1 =x2;
    y1 =y2;
  }
  return;
}

startcluster1()
{long xedge,yedge;
 int row,col,bit;
 long currentword;
 long index0,index1;
 /* flip previous cluster and adjust demons on boundary */
    for (row=0;row<NROWS;row++)
      {for (col=0;col<HWORDS;col++)
         {index0=(demonindex+2*(col+row*HWORDS));
          index1=(index0+1);
          while (index0>=NDEMONS) index0-=NDEMONS;
          while (index1>=NDEMONS) index1-=NDEMONS;
          if (0==unchecked[row][col]) /* don't bother if unchecked */
           {unchecked[row][col]=1;    /* reset unchecked */
            currentword=cluster[row][col];
             /* find cluster edges */
            xedge=currentword^((currentword<<1)
              |(1&(cluster[row][nextcol]>>(WORDSIZE-1))));
            yedge=currentword^cluster[nextrow][col];
             /* update demons and field */
            for (bit=0;bit<DEMONBITS;bit++)
             {demon[bit][index0]
               =(newdemon[bit][index0]&xedge)|((~xedge)&demon[bit][index0]);
              demon[bit][index1]
               =(newdemon[bit][index1]&yedge)|((~yedge)&demon[bit][index1]);
             }
            field[row][col]^=currentword;
             /* do any cooling or heating required */
            if (drand48()<0.05)
             {if (heater) 
               {demon[0][index0]^=RND;
                demon[0][index1]^=RND;
               }
              else if (cooler)
               {demon[0][index0]&=RND;
                demon[0][index1]&=RND;
               }
              else /* rotate the demon's first bit to help scramble */
               {demon[0][index0]=(demon[0][index0]<<1)
                   |(1&(demon[0][index0]>>(WORDSIZE-1)));
                demon[0][index1]=(demon[0][index1]<<1)
                   |(1&(demon[0][index1]>>(WORDSIZE-1)));
              }
             }
           } /* end of if checked */
            /* accumulate counts */
          currentword=demon[0][index0];
          energycount+=bitcount(currentword);
          currentword=demon[0][index1];
          energycount+=bitcount(currentword);
          volumecount+=WORDSIZE; 
         } /* end of col loop */
      } /* end of row loop */
/* clear previous cluster */
    for (col=0;col<HWORDS;col++)
      for (row=0;row<NROWS;row++)
         cluster[row][col]=0;
/* randomize starting demon location */
  demonindex=lrand48()%NDEMONS;
/* pick site to seed cluster */
  row=lrand48()%NROWS;
  col=lrand48()%HWORDS;
  clustergrowth=cluster[row][col]=1<<(lrand48()%WORDSIZE);
/* set up bounding box for cluster */
  boxleft=col+HWORDS-(1&(clustergrowth>>(WORDSIZE-1)));
  boxright=col+HWORDS+1+(clustergrowth&1);
  boxbottom=previousrow+NROWS;
  boxtop=boxbottom+3;
  /* monitor temperature */
if (volumecount>=VOLUME) /* update beta and thermometer */
 {beta=.5*log((double) (-1.+volumecount/(0.5*energycount))); 
 /* fixthermometer(); */
  energycount=0;
  volumecount=0;
 } /* end of temperature monitoring */
  /* copy lattice to window */
 XPutImage(display,window,gc,spinimage,0,0,64,132,NCOLS,NROWS); 
 return;
} /* end of startcluster */

growcluster()
/* the cluster growing routine */
{long newsites,temp;
 int row,col,rowloop,colloop,newboxbottom,newboxtop,newboxleft,newboxright;
 long currentword,temporary;
   newboxbottom=NROWS;
   newboxtop=0;
   newboxleft=HWORDS;
   newboxright=0;
   clustergrowth=0;
   direction=(-direction); /* sweep rows in alternate directions */
   for (colloop=boxleft;colloop<boxright;colloop++)
      {col=colloop;
       while (col>=HWORDS) col-=HWORDS;
       for (rowloop=boxbottom;rowloop<boxtop;rowloop++)
         {if (direction<0)
            row=(boxbottom+boxtop-rowloop-1);
          else 
            row=rowloop;
          while (row>=NROWS) row-=NROWS;
           /* make sure sadx and sady for current word are calculated */
          if (unchecked[row][col])
           check(row,col);
          if (unchecked[previousrow][col])
           check(previousrow,col);
          if (unchecked[row][previouscol])
           check(row,previouscol);
          newsites=currentword=cluster[row][col];
          temporary=sadx[row][col];
          newsites|=
            (sady[row][col]&cluster[nextrow][col]); /* up */
          newsites|=
            (sady[previousrow][col]&cluster[previousrow][col]); /* down */
          newsites|=
           (temporary&1&(cluster[row][nextcol]>>(WORDSIZE-1))); /* from right */
          newsites|=
            ((sadx[row][previouscol]&cluster[row][previouscol])
                <<(WORDSIZE-1)); /*left */
          while (newsites) /* keep working on present word */
           {newsites|=
              (temporary&(currentword<<1)); /* grow to left*/
            newsites|=
              (((temporary&currentword)>>1)&(~LEFTBIT)); /* right */
            newsites&=(~currentword);
            currentword|=newsites;
           }
              /* are new sites added? ( single "=" intentional here) */
          if (newsites=(currentword^cluster[row][col]))
             {clustergrowth=1;
              cluster[row][col]=currentword;
              rowloop-=(2*(rowloop>0)); /* back up */
              if (row<newboxbottom) 
               newboxbottom=row;
              if (row>newboxtop) 
               newboxtop=row;
              temp=col-(1&(newsites>>(WORDSIZE-1)));
              if (temp<newboxleft) newboxleft=temp;
              temp=col+(newsites&1);
              if (temp>newboxright) newboxright=temp;
             } /* end of if  new sites added */
         } /* end row loop */
       } /* end column loop */
if (clustergrowth)
 {if (newboxtop-newboxbottom<NROWS-2)
   {boxbottom=newboxbottom+NROWS-1;
    boxtop=newboxtop+NROWS+2;
   }
  else
   {boxbottom=0;
    boxtop=NROWS;
   }
  if (newboxright-newboxleft<HWORDS)
   {boxleft=newboxleft+HWORDS;
    boxright=newboxright+HWORDS+1;
   }
  else 
   {boxleft=0;
    boxright=HWORDS;
   }
   /* display cluster */
  XPutImage(display,window,gc,clusterimage,0,0,64,152+NROWS,NCOLS,NROWS); 
 }
 return;
} /* end of growcluster */

check(row,col)
int row,col;
{ /* calculate sadx and sady when unchecked */
 /* demons are sad if they cannot flip a bond */
 long right,bottom,carry0,carry1;
 int bit;
 long index0,index1;
 right=field[row][col]^((field[row][col]<<1)
      |(1&(field[row][nextcol]>>(WORDSIZE-1))));
 bottom=field[row][col]^field[nextrow][col];
 if (boundary) /* correct for antiperiodic boundaries */
  {if (row==NROWS-1)
     bottom^=(-1);
   if (col==topcol)
     right^=1;
  }  
 index0=(demonindex+2*(col+row*HWORDS)); /* start of xdemons */
 index1=(index0+1); /* start of ydemons */
 while (index0>=NDEMONS) index0-=NDEMONS;
 while (index1>=NDEMONS) index1-=NDEMONS;
 /* copy demons to newdemon */
 for (bit=0;bit<DEMONBITS;bit++)
  {newdemon[bit][index0]=demon[bit][index0];
   newdemon[bit][index1]=demon[bit][index1];
  }
 /* subtract 1 from newdemon */
 carry0=carry1=(-1);
 for (bit=0;bit<DEMONBITS;bit++)
  {newdemon[bit][index0]^=carry0;
   newdemon[bit][index1]^=carry1;
   carry0&=newdemon[bit][index0];
   carry1&=newdemon[bit][index1];
  }
 /* add two if neighbor antiparallel */
 for (bit=1;bit<DEMONBITS;bit++)
  {newdemon[bit][index0]^=right;
   newdemon[bit][index1]^=bottom;
   right&=(~newdemon[bit][index0]);
   bottom&=(~newdemon[bit][index1]);
  }
 carry0^=right;
 carry1^=bottom;
 /* make demons sad if overflow */
 sadx[row][col]=carry0;
 sady[row][col]=carry1;
 unchecked[row][col]=0;
 return;
}

setup()
{int i,j,count;
/* initialize "bits" used for bitcounts */
 for (i=0;i<256;i++)
  {j=i;
   count=j&1;
   while (j) count+=((j>>=1)&1);
   bits[i]=count; /* the number of set bits in i */
  }
 /* make checkermask have alternate bits set */
 checkermask=2;
 for (i=0;i<32;i++) checkermask|=(checkermask<<2);
 return;
}

repaint()
/* this fixes the window up whenever it is uncovered */
{int i;
 XSetPlaneMask(display,gc,AllPlanes); 
/* draw the heater buttons, taking care to erase the old active one */
/* XFillRectangle(display,heatwindow,gcr,0,
    0,HEATERWIDTH,HEATERHEIGHT);
 XDrawString(display,heatwindow,gccolor,
        0,font_height,                 "conserve E",10);
 XDrawString(display,heatwindow,gccolor,
        0,HEATERHEIGHT/3+font_height,  "   heat   ",10);
 XDrawString(display,heatwindow,gccolor,
        0,2*HEATERHEIGHT/3+font_height,"   cool   ",10);
 for (i=1;i<3;i++)
   XDrawLine(display,heatwindow,gc,
       0,i*HEATERHEIGHT/3-1,HEATERWIDTH,i*HEATERHEIGHT/3-1);
 XFillRectangle(display,heatwindow,gcxor,0,
    (HEATERHEIGHT*(heater+2*cooler))/3,HEATERWIDTH,-1+HEATERHEIGHT/3);

 XFillRectangle(display,boundwindow,gcr,0,
        0,BOUNDWIDTH,BOUNDHEIGHT);
 XDrawString(display,boundwindow,gccolor,
        0,font_height,                "  periodic  ",12);
 XDrawString(display,boundwindow,gccolor,
        0,BOUNDHEIGHT/4+font_height,  "antiperiodic",12);
 XDrawString(display,boundwindow,gccolor,
        5,2*BOUNDHEIGHT/4+font_height,"black edges",11);
 XDrawString(display,boundwindow,gccolor,
        5,3*BOUNDHEIGHT/4+font_height,"white edges",11);
 for (i=1;i<4;i++)
   XDrawLine(display,boundwindow,gc,
       0,i*BOUNDHEIGHT/4-1,BOUNDWIDTH,i*BOUNDHEIGHT/4-1);
 XFillRectangle(display,boundwindow,gcxor,0,
    (BOUNDHEIGHT*boundary)/4,BOUNDWIDTH,-1+BOUNDHEIGHT/4);

 XDrawString(display,quitbutton,gcr,
        0,font_height,"quit",4); */
/* mypause();
 myalgo();*/
/* draw thermometer*/
#define TTOP 128
#define THGT (NROWS*2)
#define TBOTTOM (TTOP+THGT)
#define TLEFT 24
#define TWDTH 4
/*
 XDrawLine(display,window,gc,TLEFT-3      ,TTOP-2   ,TLEFT+TWDTH+2,TTOP-2);
 XDrawLine(display,window,gc,TLEFT-3      ,TTOP-2   ,TLEFT-3      ,TBOTTOM+2);
 XDrawLine(display,window,gc,TLEFT-3      ,TBOTTOM+2,TLEFT+TWDTH+2,TBOTTOM+2);
   XDrawLine(display,window,gc,TLEFT+TWDTH+2,TBOTTOM+2,TLEFT+TWDTH+2,TTOP-2);
 for (i=TBOTTOM-1;i>TTOP;i-=(THGT/20))
  {XDrawLine(display,window,gc,TLEFT-7,i,TLEFT-4,i);
   XDrawLine(display,window,gc,TLEFT+TWDTH+3,i,TLEFT+TWDTH+6,i);
  }
 i=TBOTTOM-THGT*(-0.6+1.0/log(1.+sqrt(2.)));
 XDrawLine(display,window,gccolor,TLEFT-15,i,TLEFT+TWDTH+15,i); 
 XDrawString(display,window,gccolor,TLEFT+TWDTH+15,i+4,"T",1);
 XDrawString(display,window,gccolor,TLEFT+TWDTH+23,i+8,"c",1); 
  write various strings 
 sprintf(stringbuffer,"%d by %d lattice",NCOLS,NROWS);
 XDrawString(display,window,gc,140,104,stringbuffer,strlen(stringbuffer));
 XDrawString(display,window,gc,5,50,"Algorithm:",10);
 XDrawString(display,window,gc,20,104,"beta=",5);
 XDrawString(display,window,gc,64,123,"spins:",6);
 XDrawString(display,window,gc,64,147+NROWS,"changes:",8);
 XDrawString(display,window,gccolor,TLEFT+NCOLS,118,"MJC",3);   
 XPutImage(display,window,gc,spinimage,0,0,64,132,NCOLS,NROWS); 
 XPutImage(display,window,gc,clusterimage,0,0,64,152+NROWS,NCOLS,NROWS); 
 fixthermometer(); 
 fix up about screen */
 if (aboutmapped)
  {XDrawString(display,about,gc,0,2*font_height,"              XIsing    ",24);
   XDrawString(display,about,gc,0,3*font_height,"                by      ",24);
   XDrawString(display,about,gc,0,4*font_height,"          Michael Creutz",24);
   XDrawString(display,about,gc,0,5*font_height,"          creutz@bnl.gov",24);
   XDrawString(display,about,gc,0,7*font_height,     
     "   A  microcanonical simulation  of",36);
   XDrawString(display,about,gc,0,8*font_height,     
     " the two dimensional Ising model.  ",36);
   XDrawString(display,about,gc,0,10*font_height,
     " Reference for the local algorithm:",36);
   XDrawString(display,about,gc,0,11*font_height,
    " Phys. Rev. Letters 50, 1411 (1983).",36);
   XDrawString(display,about,gc,0,13*font_height,      
    " For the cluster algorithm:         ",36);
   XDrawString(display,about,gc,0,14*font_height,      
    " Phys. Rev. Letters 69, 1002 (1992).",36);
   XDrawString(display,about,gc,0,16*font_height,      
    " (click mouse to remove this window)",36);
  }
 XSetPlaneMask(display,gc,
             WhitePixel(display,screen)^BlackPixel(display,screen)); 
 return;
}

refreshdemons()
/* refresh cluster demons at beta 
To generate a bit set with probability p:
write p as a binary fraction.  Compare a
random bit string with the digits of this 
fraction and find the first place at which
they are the same.  Output the bit  
that occupies that place.  This is done here
in parallel on WORDSIZE demons at a time. */
{long morebits,nrbit;
 double bfactor,prob,temp;
 long i;
 int bit,j,pbit[16];
 bfactor=exp(-beta);
 for (bit=0;bit<DEMONBITS;bit++)
    {bfactor*=bfactor; /* gives exp(-2 beta) for the first bit, 
                                exp(-4 beta) for the next */
     prob=bfactor/(1.+bfactor); /* prob that current demon bit is 1 */
     temp=prob;
       /* calculate the binary expression for p to 16 bit accuracy */
     for (j=0;j<16;j++) 
      {temp*=2; 
       pbit[j]=(temp>1.0);
       temp-=pbit[j];
      }
     for (i=0;i<NDEMONS;i++)
      {demon[bit][i]=0;
       morebits=(-1);
       for (j=0;j<16;j++)
        {nrbit=RND;  /* new set of random bits */
         if (pbit[j]) /* comparison bit is set */
          {demon[bit][i]|=(nrbit&morebits); /* if both are equal, set demon */
           morebits&=(~nrbit); /* finished with this demon */
          }
         else /* bit of p is zero */ 
          {morebits&=nrbit; /* finished with demons with zeros in nrbit */
          }
         if (0==morebits) break; /* all bits decided */
        }
      }
    }
  return;
}

localupdate()
/* the local updating routine */
{long *up,*down; /* pointers to above and below rows */
 long temp0,temp1,reject,carry,spins,left,right,top,bottom;
 int row=0,col,rowloop;
 checkermask=(~checkermask);
 for (rowloop=0;rowloop<NROWS;rowloop++)
  {row=(row+29); /* jump ahead 29 rows */
   while (row>=NROWS) row-=NROWS;
      /* get neighbors */
   up=field[previousrow];
   down=field[nextrow];
   leftshift( field[row],work0);
   rightshift(field[row],work1);
   for (col=0;col<HWORDS;col++)
    {spins=field[row][col];
        /* determine antiparallel neighbors */
     right=work0[col]^spins;    
     left=work1[col]^spins; 
     top=up[col]^spins;
     bottom=down[col]^spins; 
     if (0==row) 
     switch(boundary)
         {case 1:
           top^=(~0);
           break;
          case 2:
           top=(~0)^spins;
           break;
          case 3:
           top=spins;
           break;
         }  
     else if (NROWS-1==row) 
        switch(boundary)
         {case 1:
           bottom^=(~0);
           break; 
          case 2:
           bottom=(~0)^spins;
           break;
          case 3:
           bottom=spins;
           break;
         }  
          /* add up top, bottom, and left antiparallel spins */
     temp0=top^bottom^left;
     temp1=(top&(bottom|left))|(bottom&left);
      /* add in right */
     carry=temp0&right; 
     temp0^=right;
     temp1^=carry;
     reject=carry&(~temp1);
      /* add  demon0 */
     carry=demon0&temp0;
     temp0^=demon0;
     temp1^=carry;
     reject^=(carry&(~temp1));
      /* add in demon1 */
     temp1^=demon1;
     reject^=(demon1&(~temp1));
      /*  subtract two */
     temp1=~temp1; 
     reject^=temp1;
      /* only change one color sites */
     reject|=checkermask;
      /* make changes */
     demon0=(demon0&reject)|(temp0&(~reject));
     demon1=(demon1&reject)|(temp1&(~reject));
     field[row][col]=(spins^(~reject)); 
     cluster[row][col]&=checkermask; 
     cluster[row][col]|=(~reject); 
     energycount+=bitcount(demon0&(~checkermask));
    }
           /* heat or cool */
   if (heater) 
       {if (drand48()<.04) demon0=RND;
       }
   else if (cooler) 
       {if (drand48()<.04) demon0=demon1=0;
       }
  }
 iteration=(iteration+1);
   /* monitor temperature every 10 "half" iterations */
 if (iteration==10) /* update beta and thermometer */
  {iteration=0;
   beta=.25*log((double) (-1.+10*VOLUME/(2.*energycount)));
 /*  fixthermometer(); */
   energycount=0;
  }
  /* copy lattice to window */
 XPutImage(display,window,gc,spinimage,0,0,64,132,NCOLS,NROWS); 
/* XPutImage(display,window,gc,clusterimage,0,0,64,152+NROWS,NCOLS,NROWS); */
 return;
}
 
fixthermometer()
{int temp;
 sprintf(stringbuffer,"%.3f",beta);
 temp=THGT*(-0.6+0.5/beta); 
 if (temp<0) temp=0;
 if  (temp>THGT) temp=THGT;
 XDrawImageString(display,window,gccolor,64,104
   ,stringbuffer,strlen(stringbuffer));
 XFillRectangle(display,window,gccolor,TLEFT,TBOTTOM-temp,TWDTH,temp);
 XFillRectangle(display,window,gcr,TLEFT,TTOP,TWDTH,THGT-temp);
 return;
}

rightshift(source,destination)
/* shifts a row of cells to the right by one bit */
long source[HWORDS],destination[HWORDS];
{int col;
 long carry=0;
 for (col=bottomcol;col!=(topcol+nshift);col+=nshift)
  {destination[col]=(~LEFTBIT)&(source[col]>>1);
   if (carry) destination[col]|=LEFTBIT;
   carry=source[col]&1;
  }
 if (carry) destination[bottomcol]|=LEFTBIT;
 switch(boundary)
   {case 1:
     destination[bottomcol]^=LEFTBIT;
     break;
    case 2:
     destination[bottomcol]|=LEFTBIT;
     break;
    case 3:
     destination[bottomcol]&=(~LEFTBIT);
     break;
   }
 return;
}

leftshift(source,destination)
/* shifts a row of cells one bit to the left */
long source[HWORDS],destination[HWORDS];
{int col;
 long carry=0;
 for (col=topcol;col!=(bottomcol-nshift);col-=nshift)
  {destination[col]=(source[col]<<1);
   if (carry) destination[col]|=1;
   carry=source[col]>>(WORDSIZE-1);
  }
 if (carry) destination[topcol]|=1;
 switch(boundary)
   {case 1:
     destination[topcol]^=1;
     break;
   case 2:
     destination[topcol]|=1;
     break;
   case 3:
     destination[topcol]&=(~1);
     break;
   }
 return;
} 

/* pause the updating */
mypause()
{if (0==paused) 
   XDrawImageString(display,pausebutton,gcr,
        0,font_height,"pause",5);
 else
   XDrawImageString(display,pausebutton,gcr,
        0,font_height," run ",5);
 return;
}

myalgo()
{if (0==algorithm) 
   XDrawImageString(display,algobutton,gcr,
        0,font_height," local ",7);
 else
   XDrawImageString(display,algobutton,gcr,
        0,font_height,"cluster",7);
 return;
}

setalg()
/* take action if algorithm button hit*/
{int row,col;
 energycount=volumecount=iteration=0;
 if (1==algorithm) /* reset cluster stuff */
  {if (boundary>1)
     {setboundary(0);
     }
   refreshdemons(); 
   clustergrowth=0;
   for (col=0;col<HWORDS;col++)
   for (row=0;row<NROWS;row++)
    {cluster[row][col]=0;
     unchecked[row][col]=1;
    }
  }
 return;
}

setheater(value)
/* acts on heater buttons */
int value;
{value&=3;
/* XFillRectangle(display,heatwindow,gcxor,0,
    (HEATERHEIGHT*(heater+2*cooler))/3,HEATERWIDTH,-1+HEATERHEIGHT/3);
 heater=1&value;
 cooler=value>>1;
 XFillRectangle(display,heatwindow,gcxor,0,
    (HEATERHEIGHT*(heater+2*cooler))/3,HEATERWIDTH,-1+HEATERHEIGHT/3);
*/
 return;
}

setboundary(value)
/* take action if boundary buttons hit */
int value;
{XFillRectangle(display,boundwindow,gcxor,0,
    (BOUNDHEIGHT*boundary)/4,BOUNDWIDTH,-1+BOUNDHEIGHT/4);
/* don't allow fixed boundaries with cluster algorithm */
 if ((algorithm==1)&&(value>1))
   {fprintf(stderr,"fixed boundaries not implemented for cluster algorithm\n");
   }
 else
  boundary=value;
 XFillRectangle(display,boundwindow,gcxor,0,
    (BOUNDHEIGHT*boundary)/4,BOUNDWIDTH,-1+BOUNDHEIGHT/4);
 return;
}

openwindow(argc,argv)
/* a lot of this is taken from the basicwin program in the
Xlib Programming Manual */
int argc;
char **argv;
{char *window_name="Microcanonical Ising Simulation";
 char *icon_name="Ising";
 Pixmap icon_pixmap;
 char *display_name=NULL;
 long event_mask;
 XColor xcolor,colorcell;
 Colormap cmap;
 int depth,darkcolor,lightcolor;

#define icon_bitmap_width 16
#define icon_bitmap_height 16
 static char icon_bitmap_bits[] = {
   0x1f, 0xf8, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0x88, 0x1f, 0xf8,
   0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0x1f, 0xf8, 0xff, 0xff,
   0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};

/* open up the display */
 if ((display=XOpenDisplay(display_name))==NULL)
  {fprintf(stderr,"%s: cannot connect to X server %s\n",
    progname,XDisplayName(display_name));
   exit(-1);
  }
 screen=DefaultScreen(display);

 depth=DefaultDepth(display,screen);
 cmap=DefaultColormap(display,screen);

 darkcolor=BlackPixel(display,screen);
 lightcolor=WhitePixel(display,screen);
 if (depth>1) /* color? This is not the right way to do it, but .... */
  {if (XAllocNamedColor(display,cmap,"firebrick",&colorcell,&xcolor))
              darkcolor=colorcell.pixel;
   if (XAllocNamedColor(display,cmap,"wheat",&colorcell,&xcolor))
              lightcolor=colorcell.pixel;
  }

/* make the main window */
 windowwidth=NCOLS+75;
 windowheight=2*NROWS+164;
 windowwidth=windowheight=800;
 window=XCreateSimpleWindow(display,RootWindow(display,screen),
   0,0,windowwidth,windowheight,4,BlackPixel(display,screen),
   lightcolor);

/* make the icon */
 icon_pixmap=XCreateBitmapFromData(display,window,
   icon_bitmap_bits,icon_bitmap_width,
   icon_bitmap_height);

 size_hints.flags=PPosition | PSize | PMinSize;
 size_hints.min_width=windowwidth;
 size_hints.min_height=windowheight;
#ifdef X11R3
 size_hints.x=x;
 size_hints.y=y;
 size_hints.width=windowwidth;
 size_hints.height=windowheight;
 XSetStandardProperties(display,window,window_name,icon_name,
    icon_pixmap,argv,argc,&size_hints);
#else
 {XWMHints wm_hints;
  XClassHint class_hints;
  XTextProperty windowName, iconName;
  if (XStringListToTextProperty(&window_name,1,&windowName)==0)
   {fprintf(stderr,"%s: structure allocation for windowName failed.\n"
      ,progname);
    exit(-1);
   }
  if (XStringListToTextProperty(&icon_name,1,&iconName)==0)
   {fprintf(stderr,"%s: structure allocation for iconName failed.\n"
       ,progname);
    exit(-1);
   }
  wm_hints.initial_state=NormalState;
  wm_hints.input=True;
  wm_hints.icon_pixmap=icon_pixmap;
  wm_hints.flags=StateHint|IconPixmapHint|InputHint;
  class_hints.res_name=progname;
  class_hints.res_class="Basicwin";
  XSetWMProperties(display,window,&windowName,&iconName,
       argv,argc,&size_hints,&wm_hints,&class_hints);
 }
#endif

/* make the buttons */
 quitbutton=XCreateSimpleWindow(display,window,
    6,0,40,20,2,BlackPixel(display,screen),
    darkcolor);
 pausebutton=XCreateSimpleWindow(display,window,
    56,0,48,20,2,BlackPixel(display,screen),
    darkcolor);
 algobutton=XCreateSimpleWindow(display,window,
    8,55,64,20,2,BlackPixel(display,screen),
    darkcolor);
 heatwindow=XCreateSimpleWindow(display,window,
    windowwidth-HEATERWIDTH-BOUNDWIDTH-12,0,
    HEATERWIDTH,HEATERHEIGHT-1,2,BlackPixel(display,screen),
    lightcolor);
 boundwindow=XCreateSimpleWindow(display,window,
             windowwidth-BOUNDWIDTH-4,0,BOUNDWIDTH,BOUNDHEIGHT-1,2,
             BlackPixel(display,screen),lightcolor);
 about=XCreateSimpleWindow(display,window,
    0,3*windowheight/5,windowwidth-4,2*windowheight/5-4,2,
    BlackPixel(display,screen),WhitePixel(display,screen));

/* pick the events to look for */
 event_mask=ExposureMask|ButtonPressMask|StructureNotifyMask;
 XSelectInput(display,window,event_mask);
 XSelectInput(display,about,event_mask);
 event_mask=ButtonPressMask; 
/* note that with this simple mask if one just covers a button
it will not get redrawn.  I wonder if anyone will notice?  If I put
the exposuremask in here, things flash irritatingly on being uncovered. */
 XSelectInput(display,quitbutton,event_mask);
 XSelectInput(display,pausebutton,event_mask);
 XSelectInput(display,algobutton,event_mask);
 XSelectInput(display,heatwindow,event_mask);
 XSelectInput(display,boundwindow,event_mask);

/* pick font: 9x15 is supposed to almost always be there */
 if ((font=XLoadQueryFont(display,"9x15"))==NULL)
  {fprintf(stderr,"%s: Cannot open 9x15 font\n",progname);
   exit(-1);
  }
 font_height=font->ascent+font->descent;

/* make graphics contexts: 
      gc for black on white
      gccolor for background and buttons 
      gcr for reverse video 
      gcxor for highlighting */
   
 gc=XCreateGC(display,window,0,NULL);
 XSetFont(display,gc,font->fid);
 XSetForeground(display,gc,BlackPixel(display,screen));
 XSetBackground(display,gc,WhitePixel(display,screen)); 

 gcr=XCreateGC(display,window,0,NULL); 
 XSetFont(display,gcr,font->fid);
 XSetForeground(display,gcr,lightcolor);
 XSetBackground(display,gcr,darkcolor); 

 gccolor=XCreateGC(display,window,0,NULL); 
 XSetFont(display,gccolor,font->fid);
 XSetForeground(display,gccolor,darkcolor);
 XSetBackground(display,gccolor,lightcolor); 

 gcxor=XCreateGC(display,window,0,NULL);
 XSetFont(display,gcxor,font->fid);
 XSetForeground(display,gcxor,
             darkcolor^lightcolor);
 XSetFunction(display,gcxor,GXxor);

/* show the window and buttons */
 XMapWindow(display,window);
/* XMapWindow(display,quitbutton);
 XMapWindow(display,pausebutton);
 XMapWindow(display,algobutton);
 XMapWindow(display,heatwindow);
 XMapWindow(display,boundwindow);
*/
/* make image structures */
 spinimage=XCreateImage(display,(Visual *) &window,1,XYBitmap,0,
    (char*) field,NCOLS,NROWS,32,0);
 clusterimage=XCreateImage(display,(Visual *) &window,1,XYBitmap,0,
    (char*) cluster,NCOLS,NROWS,32,0);
 (*spinimage).bitmap_unit=WORDSIZE;
 (*clusterimage).bitmap_unit=WORDSIZE; 

/* test for byte order on client */
 field[0][0]=1;
 if ( * (char*) field)
  {(*spinimage).byte_order=LSBFirst;
   (*clusterimage).byte_order=LSBFirst;
     /* printf("byte_order=LSBFirst\n"); */
  }
 else
  {(*spinimage).byte_order=MSBFirst;
   (*clusterimage).byte_order=MSBFirst;
    /* printf("byte_order=MSBFirst\n"); */
  }
  /* find out in what order bits are put on screen */
 field[0][0]=0xff;
 if (XGetPixel(spinimage,0,0)) /* bits are backward */
    {nshift=(-1);
     topcol=0;
     bottomcol=HWORDS-1;
     (*spinimage).bitmap_bit_order=LSBFirst;
     (*clusterimage).bitmap_bit_order=LSBFirst;
     /* printf("bitmap_bit_order=LSBFirst\n"); */
    }
 else 
    {nshift=1;
     topcol=HWORDS-1;
     bottomcol=0;
     (*spinimage).bitmap_bit_order=MSBFirst;
     (*clusterimage).bitmap_bit_order=MSBFirst;
    }
 return;
}

