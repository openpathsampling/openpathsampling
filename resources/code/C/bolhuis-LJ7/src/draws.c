#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "glut.h"
#include "path.h"

typedef struct atom {

  float radius, length, theta, phi;

  vector r, u;
  int no;
} Atom;

Slice *currentslice;
int Npart;         						/* Number of particles */
int Precision = 4; 						/* Number of sides per sphere, etc */
int Quit = FALSE;  						/* Quit flag */
int Run = FALSE;   						/* Run flag */
int Drawall = 3;   						/* Draw flag*/

Atom *pos;           					/* Storage of positions, etc */
vector Boxl1, Boxl2; 					/* Size of the periodic box */
vector *Circle = NULL;

float Far = 200.0; 						/* The far clipping plane */
float Boxdistance = 1.25;
float Distance = -10;   				/* Distance, Azimuth, Inclination and Twist */
GLfloat xRotation = 90, 				/* define the viewing position */
    yRotation = 0, zRotation = 0;
char string[] = "A";

struct glstattype {

  int Animate, BackwardAnimation, OneUp, OneDown;

} glstat;

static float MaterialRed[] = { 0.7, 0.0, 0.0, 1.0 };
static float MaterialGreen[] = { 0.1, 0.5, 0.2, 1.0 };
static float MaterialBlue[] = { 0.0, 0.0, 0.7, 1.0 };
static float MaterialSome[] = { 0.0, 0.7, 0.0, 1.0 };
static float MaterialYellow[] = { 0.7, 0.7, 0.0, 1.0 };
static float MaterialCyan[] = { 0., 0.7, 0.7, 1.0 };
static float MaterialMagenta[] = { 0.7, 0., 0.7, 1.0 };

static float front_shininess[] = { 60.0 };
static float front_specular[] = { 0.7, 0.7, 0.7, 1.0 };
static float ambient[] = { 0.0, 0.0, 0.0, 1.0 };
static float diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
static float position0[] = { 1.0, 1.0, 1.0, 0.0 };
static float position1[] = { -1.0, -1.0, 1.0, 0.0 };
static float lmodel_ambient[] = { 0.5, 0.5, 0.5, 1.0 };
static float lmodel_twoside[] = { GL_TRUE };

static float vert[17][3] = { { -1.0, -1.0, -1.0 },
                             { -1.0, -1.0, 1.0 },
                             { 1.0, -1.0, 1.0 },
                             { 1.0, -1.0, -1.0 },
                             { -1.0, -1.0, -1.0 },
                             { -1.0, 1.0, -1.0 },
                             { 1.0, 1.0, -1.0 },
                             { 1.0, -1.0, -1.0 },
                             { 1.0, -1.0, 1.0 },
                             { 1.0, 1.0, 1.0 },
                             { 1.0, 1.0, 1.0 },
                             { 1.0, 1.0, -1.0 },
                             { 1.0, 1.0, 1.0 },
                             { -1.0, 1.0, 1.0 },
                             { -1.0, 1.0, -1.0 },
                             { -1.0, 1.0, 1.0 },
                             { -1.0, -1.0, 1.0 }, };

static void DrawBitmapString(void *font, const char *string) {
  int i;

  for (i = 0; string[i]; i++)
    glutBitmapCharacter(font, string[i]);
}

static void DrawStrokeString(void *font, const char *string) {
  int i;

  for (i = 0; string[i]; i++)
    glutStrokeCharacter(font, string[i]);
}

void convert_coordinates(Slice *psl) {
  int ibox, i, j, ii = 0;
  double dr2, l, skin, a, rui, lo2, rm2, rx2, ryz, d;
  vector r, r0, u, dr, ri, rc[4];

  Npart = sys.npart;
  Boxl1 = sys.boxl;
  Boxl2 = sys.boxl;

  for (i = 0; i < Npart; i++) {
    r = psl->pts[i].r;
    u = nulvec;
    pos[i].r = r;

    pos[i].u = u;
    pos[i].no = 1;
    pos[i].theta = 180.0 / PI * acos(u.z);
    if (u.x > 0) {
      pos[i].phi = 180.0 / PI * atan(u.y / u.x);
    } else
      pos[i].phi = 90;
    pos[i].length = 0;
    pos[i].radius = 0.5;
  }
}
void display_sph1(int ind) {

  glPushMatrix();
  glTranslatef(pos[ind].r.x, pos[ind].r.y, pos[ind].r.z);
  glScalef(pos[ind].radius, pos[ind].radius, pos[ind].radius);
  switch (ind) {
  case 0:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialCyan);
    break;
  case 1:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialRed);
    break;
  case 2:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialYellow);
    break;
  case 3:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialGreen);
    break;
  case 4:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialBlue);
    break;
  case 5:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialMagenta);
    break;
  case 6:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialSome);
    break;
  case 7:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialCyan);
    break;
  case 8:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialRed);
    break;
  case 9:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialYellow);
    break;
  case 10:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialGreen);
    break;
  case 11:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialBlue);
    break;
  case 12:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialMagenta);
    break;
  case 13:
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialSome);
    break;
  }
  glutSolidSphere(1, Precision, Precision);

  glPopMatrix();
}

void create_circle(void) {
  vector *ptch;
  float alpha, delta;
  int i, j;

  delta = 2.0 * PI / (16. * Precision);
  if (Circle)
    free(Circle);
  Circle = (vector *)malloc((16 * Precision + 1) * sizeof(vector));
  ptch = Circle;

  alpha = 0.0;
  for (i = 0; i < 16 * Precision; i++) {
    ptch[i].x = 0;
    ptch[i].y = sin(alpha);
    ptch[i].z = cos(alpha);
    alpha += delta;
  }
}

void convert(vector *p, float f[3]) {
  f[0] = p->x;
  f[1] = p->y;
  f[2] = p->z;
  return;
}

void lighted_plate(void) {

  vector *ptch;
  int i, j, k;
  float f[3];

  ptch = Circle;
  {
    glBegin(GL_POLYGON);
    for (k = 0; k < Precision * 16; k++) {
      glVertex3f(ptch->x, ptch->y, ptch->z);
      glNormal3f(ptch->x, ptch->y, ptch->z);
      ptch++;
    }
    glEnd();
  }

  ptch--;
  {
    glBegin(GL_POLYGON);

    for (k = 0; k < Precision * 16; k++) {
      glVertex3f(ptch->x, ptch->y, ptch->z);
      glNormal3f(ptch->x, ptch->y, ptch->z);
      ptch--;
    }
    glEnd();
  }
}

void display_plates(int ind) {

  glPushMatrix();
  glTranslatef(pos[ind].r.x, pos[ind].r.y, pos[ind].r.z);
  glScalef(1, pos[ind].radius, pos[ind].radius);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialGreen);
  lighted_plate();
  glPopMatrix();
}

void draw_box(float *vec, float boxdis) {
  int i;

  glPushMatrix();
  glScalef(0.5 * vec[0], 0.5 * vec[1], 0.5 * vec[2]);
  glTranslatef(0, boxdis, 0);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialRed);
  glColor3f(1.0, 0, 0);
  glBegin(GL_LINE_STRIP);
  for (i = 0; i < 17; i++) {
    glVertex3fv(vert[i]);
  }

  glEnd();
  glPopMatrix();
}

static void Reshape(int width, int height) {
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(-1.0, 1.0, -1.0, 1.0, 5.0, 100.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(0.0, 0.0, -15.0);
}

void Display(void) {
  float vec[3];
  int i;

  glLoadIdentity();
  glTranslatef(0.0, 0.0, Distance);
  glRotatef(xRotation, 1, 0, 0);
  glRotatef(yRotation, 0, 1, 0);
  glRotatef(zRotation, 0, 0, 1);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  vec[0] = Boxl1.x;
  vec[1] = Boxl1.y;
  vec[2] = Boxl1.z;
  draw_box(vec, 0);

  for (i = 0; i < Npart; i++) {
    if (Drawall == 3)
      display_sph1(i);
    if ((Drawall == 2) && (pos[i].no >= 2))
      display_sph1(i);
    if ((Drawall <= 1) && (pos[i].no == 3))
      display_sph1(i);
  }

  glFlush();

  glutSwapBuffers();
}

static void CursorKeys(int key, int x, int y) {

  switch (key) {
  case GLUT_KEY_LEFT:
    zRotation += 5;
    break;
  case GLUT_KEY_RIGHT:
    zRotation -= 5;
    break;
  case GLUT_KEY_UP:
    xRotation += 5;
    break;
  case GLUT_KEY_DOWN:
    xRotation -= 5;
    break;
  default:
    return;
  }
  glutPostRedisplay();
}
static void Key(unsigned char key, int x, int y) {

  switch (key) {
  case 27:
    exit(1);

  case ',':
    yRotation += 5;
    break;
  case '.':
    yRotation -= 5;
    break;
  case 'x':
    Distance += 2;
    break;
  case 'z':
    Distance -= 2;
    break;
  case '1':
    Precision = 2;
    create_circle();
    break;
  case '2':
    Precision = 4;
    create_circle();
    break;
  case '3':
    Precision = 8;
    create_circle();
    break;
  case '4':
    Precision = 16;
    create_circle();
    break;
  case '5':
    Precision = 32;
    create_circle();
    break;
  case 'l':
    Drawall++;
    if (Drawall > 3)
      Drawall = 0;
    break;
  case 'm':
    glstat.Animate = !glstat.Animate;
    break;
  case 'b':
    glstat.BackwardAnimation = !glstat.BackwardAnimation;
    break;
  case '=':
    glstat.OneUp = 1;
    break;
  case '-':
    glstat.OneDown = 1;
    break;
  default:
    return;
  }
  glutPostRedisplay();
}

static void Init(void) {
  static GLfloat mat[4] = { 0.8, 0.8, 0.0, 1.0 };
  static GLfloat pos[4] = { -1.0, 1.0, 1.0, 0.0 };

  /* setup lighting, etc */
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, position0);
  glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT1, GL_POSITION, position1);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);

  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, front_shininess);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, front_specular);

  glEnable(GL_CULL_FACE);

  glDisable(GL_RESCALE_NORMAL_EXT);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);
  glHint(GL_FOG_HINT, GL_FASTEST);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);

  create_circle();
}

#define UNSCALED 1
#define NORMALIZE 2
#define RESCALE 3
#define QUIT 4

static void ModeMenu(int entry) {
  if (entry == UNSCALED) {
    glDisable(GL_RESCALE_NORMAL_EXT);
    glDisable(GL_NORMALIZE);
  } else if (entry == NORMALIZE) {
    glEnable(GL_NORMALIZE);
    glDisable(GL_RESCALE_NORMAL_EXT);
  } else if (entry == RESCALE) {
    glDisable(GL_NORMALIZE);
    glEnable(GL_RESCALE_NORMAL_EXT);
  } else if (entry == QUIT) {
    exit(0);
  }
  glutPostRedisplay();
}

static void Idle() {
  if (sys.sim_type == 1) {
    if (glstat.Animate) {
      if (glstat.BackwardAnimation) {
        glstat.OneDown = glstat.Animate;
      } else
        glstat.OneUp = glstat.Animate;
    }

    if (glstat.OneDown) {
		currentslice--;
      if (currentslice < slice) {
        currentslice = &slice[path.nslices - 1];
      }
      glstat.OneDown = 0;
      convert_coordinates(currentslice);
      glutPostRedisplay();
    }
    if (glstat.OneUp) {
		currentslice += 10;
      if (currentslice > &slice[path.nslices - 1]) {
        currentslice = &slice[0];
        mainloop_for_graphics();
      }
      glstat.OneUp = 0;
      convert_coordinates(currentslice);
      glutPostRedisplay();
    }

    if ((glstat.OneDown == 0) && (glstat.OneUp == 0)) {
    }
  } else {
    if (glstat.Animate) {
      glstat.OneUp = glstat.Animate;
    }
    if (glstat.OneUp) {
      if (sys.sim_type == 3) {
        mainloop_for_graphics();
        currentslice = &slice[0];
      }
      convert_coordinates(currentslice);
      glstat.OneUp = 0;
      glutPostRedisplay();
    }
  }

  return;
}

void Init_Graphics(int argc, char *argv[], Slice *psl) {

  printf("In graphics\n");

  glutInit(&argc, argv);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(640, 800);

  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);

  glutCreateWindow(argv[0]);
  glClearDepth(1.0);
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glColor3f(1.0, 1.0, 1.0);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glFlush();
  glutSwapBuffers();

  Init();

  Npart = sys.npart;
  pos = (Atom *)calloc(Npart, sizeof(Atom));

  convert_coordinates(psl);
  currentslice = psl;

  glutIdleFunc(Idle);
  glutReshapeFunc(Reshape);
  glutSpecialFunc(CursorKeys);
  glutKeyboardFunc(Key);
  glutDisplayFunc(Display);

  glutCreateMenu(ModeMenu);
  glutAddMenuEntry("Unscaled", UNSCALED);
  glutAddMenuEntry("Normalize", NORMALIZE);
  glutAddMenuEntry("Rescale EXT", RESCALE);
  glutAddMenuEntry("Quit", QUIT);
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  glutMainLoop();
  return;
}
