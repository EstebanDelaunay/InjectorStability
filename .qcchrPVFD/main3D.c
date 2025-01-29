#line 1 "main3D-cpp.c"
#line 1 "<built-in>"
#line 1 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 1 "<command-line>"
#line 1 "main3D-cpp.c"
#if _XOPEN_SOURCE < 700
  #undef _XOPEN_SOURCE
  #define _XOPEN_SOURCE 700
#endif
#if _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif



#line 1 "common.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/common.h"





#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define sign2(x) ((x) > 0 ? 1 : (x) < 0 ? -1 : 0)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))

#define unmap(x,y)

#define trash(x)


#line 1 "grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/stdlib.h"
#include <stdlib.h>
#line 2 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 3 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/stddef.h"
#include <stddef.h>
#line 4 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/stdbool.h"
#include <stdbool.h>
#line 5 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/stdarg.h"
#include <stdarg.h>
#line 6 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/string.h"
#include <string.h>
#line 7 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/float.h"
#include <float.h>
#line 8 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/limits.h"
#include <limits.h>
#line 9 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/math.h"
#include <math.h>
#line 10 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/time.h"
#include <time.h>
#line 11 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/sys/time.h"
#include <sys/time.h>
#line 12 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/sys/resource.h"
#include <sys/resource.h>
#line 13 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 22 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
# define OMP(x)

# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe
#line 48 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : 1e30)

#define systderr stderr
#define systdout stdout

FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 65

# define system(command) (pid() == 0 ? system(command) : 0)
#line 76 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (ferr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}







#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup




# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)






#line 1 "grid/array.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/array.h"


typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void * array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
  return (void *)(((char *)a->p) + a->len - size);
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
  return p;
}
#line 108 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#line 351 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
# define tracing(...)
# define end_tracing(...)
#line 367 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  if (!(!in_prof)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/config.h", 370, "!in_prof"); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 372

#define prof_stop()\
  if (!(in_prof)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/config.h", 374, "in_prof"); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 377






     
int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{tracing("mpi_all_reduce0","/home/esteban/ProgramFile/basilisk/src/grid/config.h",384);
  { int _ret= MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm);end_tracing("mpi_all_reduce0","/home/esteban/ProgramFile/basilisk/src/grid/config.h",387);return _ret;}
end_tracing("mpi_all_reduce0","/home/esteban/ProgramFile/basilisk/src/grid/config.h",388);}
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 396

#define mpi_all_reduce_array(v,type,op,elem) {\
  prof_start ("mpi_all_reduce");\
  type * global = malloc ((elem)*sizeof(type)), * tmp = malloc ((elem)*sizeof(type));\
  for (int i = 0; i < elem; i++)\
    tmp[i] = (v)[i];\
  MPI_Datatype datatype;\
  if (!strcmp(#type, "double")) datatype = MPI_DOUBLE;\
  else if (!strcmp(#type, "int")) datatype = MPI_INT;\
  else if (!strcmp(#type, "long")) datatype = MPI_LONG;\
  else if (!strcmp(#type, "bool")) datatype = MPI_C_BOOL;\
  else if (!strcmp(#type, "unsigned char")) datatype = MPI_UNSIGNED_CHAR;\
  else {\
    fprintf (stderr, "unknown reduction type '%s'\n", #type);\
    fflush (stderr);\
    abort();\
  }\
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);\
  for (int i = 0; i < elem; i++)\
    (v)[i] = global[i];\
  free (global), free (tmp);\
  prof_stop();\
}\

#line 419




#define QFILE FILE

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#line 500 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
  }
}
#line 513 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#define OMP_PARALLEL() OMP(omp parallel)

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define _index(a,m) (a.i)
#define val(a,k,l,m) data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;
#line 530 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"
#if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA && !_GPU
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  if (!(sizeof (int64_t) == sizeof (double))) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/config.h", 540, "sizeof (int64_t) == sizeof (double)");
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif



static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}

#line 1 "grid/../ast/symbols.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/../ast/symbols.h"

const char * symbol_name (int sym);
int token_symbol (int token);
enum yysymbol_kind_t
{
  sym_YYEMPTY = -2,
  sym_YYEOF = 0,
  sym_YYerror = 1,
  sym_YYUNDEF = 2,
  sym_IDENTIFIER = 3,
  sym_I_CONSTANT = 4,
  sym_F_CONSTANT = 5,
  sym_STRING_LITERAL = 6,
  sym_FUNC_NAME = 7,
  sym_SIZEOF = 8,
  sym_PTR_OP = 9,
  sym_INC_OP = 10,
  sym_DEC_OP = 11,
  sym_LEFT_OP = 12,
  sym_RIGHT_OP = 13,
  sym_LE_OP = 14,
  sym_GE_OP = 15,
  sym_EQ_OP = 16,
  sym_NE_OP = 17,
  sym_AND_OP = 18,
  sym_OR_OP = 19,
  sym_MUL_ASSIGN = 20,
  sym_DIV_ASSIGN = 21,
  sym_MOD_ASSIGN = 22,
  sym_ADD_ASSIGN = 23,
  sym_SUB_ASSIGN = 24,
  sym_LEFT_ASSIGN = 25,
  sym_RIGHT_ASSIGN = 26,
  sym_AND_ASSIGN = 27,
  sym_XOR_ASSIGN = 28,
  sym_OR_ASSIGN = 29,
  sym_TYPEDEF_NAME = 30,
  sym_ENUMERATION_CONSTANT = 31,
  sym_TYPEDEF = 32,
  sym_EXTERN = 33,
  sym_STATIC = 34,
  sym_AUTO = 35,
  sym_REGISTER = 36,
  sym_INLINE = 37,
  sym_CONST = 38,
  sym_RESTRICT = 39,
  sym_VOLATILE = 40,
  sym_BOOL = 41,
  sym_CHAR = 42,
  sym_SHORT = 43,
  sym_INT = 44,
  sym_LONG = 45,
  sym_SIGNED = 46,
  sym_UNSIGNED = 47,
  sym_FLOAT = 48,
  sym_DOUBLE = 49,
  sym_VOID = 50,
  sym_COMPLEX = 51,
  sym_IMAGINARY = 52,
  sym_STRUCT = 53,
  sym_UNION = 54,
  sym_ENUM = 55,
  sym_ELLIPSIS = 56,
  sym_CASE = 57,
  sym_DEFAULT = 58,
  sym_IF = 59,
  sym_ELSE = 60,
  sym_SWITCH = 61,
  sym_WHILE = 62,
  sym_DO = 63,
  sym_FOR = 64,
  sym_GOTO = 65,
  sym_CONTINUE = 66,
  sym_BREAK = 67,
  sym_RETURN = 68,
  sym_ALIGNAS = 69,
  sym_ALIGNOF = 70,
  sym_ATOMIC = 71,
  sym_GENERIC = 72,
  sym_NORETURN = 73,
  sym_STATIC_ASSERT = 74,
  sym_THREAD_LOCAL = 75,
  sym_MAYBECONST = 76,
  sym_NEW_FIELD = 77,
  sym_TRACE = 78,
  sym_FOREACH = 79,
  sym_FOREACH_INNER = 80,
  sym_FOREACH_DIMENSION = 81,
  sym_REDUCTION = 82,
  sym_83_ = 83,
  sym_84_ = 84,
  sym_85_ = 85,
  sym_86_ = 86,
  sym_87_ = 87,
  sym_88_ = 88,
  sym_89_ = 89,
  sym_90_ = 90,
  sym_91_ = 91,
  sym_92_ = 92,
  sym_93_ = 93,
  sym_94_ = 94,
  sym_95_ = 95,
  sym_96_ = 96,
  sym_97_ = 97,
  sym_98_ = 98,
  sym_99_ = 99,
  sym_100_ = 100,
  sym_101_ = 101,
  sym_102_ = 102,
  sym_103_ = 103,
  sym_104_ = 104,
  sym_105_ = 105,
  sym_106_ = 106,
  sym_YYACCEPT = 107,
  sym_translation_unit = 108,
  sym_primary_expression = 109,
  sym_expression_error = 110,
  sym_constant = 111,
  sym_enumeration_constant = 112,
  sym_string = 113,
  sym_generic_selection = 114,
  sym_generic_assoc_list = 115,
  sym_generic_association = 116,
  sym_postfix_expression = 117,
  sym_postfix_initializer = 118,
  sym_array_access = 119,
  sym_function_call = 120,
  sym_member_identifier = 121,
  sym_argument_expression_list = 122,
  sym_argument_expression_list_item = 123,
  sym_unary_expression = 124,
  sym_unary_operator = 125,
  sym_cast_expression = 126,
  sym_multiplicative_expression = 127,
  sym_additive_expression = 128,
  sym_shift_expression = 129,
  sym_relational_expression = 130,
  sym_equality_expression = 131,
  sym_and_expression = 132,
  sym_exclusive_or_expression = 133,
  sym_inclusive_or_expression = 134,
  sym_logical_and_expression = 135,
  sym_logical_or_expression = 136,
  sym_conditional_expression = 137,
  sym_assignment_expression = 138,
  sym_assignment_operator = 139,
  sym_expression = 140,
  sym_constant_expression = 141,
  sym_declaration = 142,
  sym_declaration_specifiers = 143,
  sym_init_declarator_list = 144,
  sym_init_declarator = 145,
  sym_storage_class_specifier = 146,
  sym_type_specifier = 147,
  sym_types = 148,
  sym_struct_or_union_specifier = 149,
  sym_struct_or_union = 150,
  sym_struct_declaration_list = 151,
  sym_struct_declaration = 152,
  sym_specifier_qualifier_list = 153,
  sym_struct_declarator_list = 154,
  sym_struct_declarator = 155,
  sym_enum_specifier = 156,
  sym_enumerator_list = 157,
  sym_enumerator = 158,
  sym_atomic_type_specifier = 159,
  sym_type_qualifier = 160,
  sym_function_specifier = 161,
  sym_alignment_specifier = 162,
  sym_declarator = 163,
  sym_direct_declarator = 164,
  sym_generic_identifier = 165,
  sym_pointer = 166,
  sym_type_qualifier_list = 167,
  sym_parameter_type_list = 168,
  sym_parameter_list = 169,
  sym_parameter_declaration = 170,
  sym_identifier_list = 171,
  sym_type_name = 172,
  sym_abstract_declarator = 173,
  sym_direct_abstract_declarator = 174,
  sym_type_not_specified = 175,
  sym_initializer = 176,
  sym_initializer_list = 177,
  sym_designation = 178,
  sym_designator_list = 179,
  sym_designator = 180,
  sym_static_assert_declaration = 181,
  sym_statement = 182,
  sym_labeled_statement = 183,
  sym_compound_statement = 184,
  sym_185_1 = 185,
  sym_block_item_list = 186,
  sym_block_item = 187,
  sym_expression_statement = 188,
  sym_selection_statement = 189,
  sym_for_scope = 190,
  sym_iteration_statement = 191,
  sym_for_declaration_statement = 192,
  sym_jump_statement = 193,
  sym_external_declaration = 194,
  sym_function_declaration = 195,
  sym_function_definition = 196,
  sym_declaration_list = 197,
  sym_basilisk_statements = 198,
  sym_macro_statement = 199,
  sym_foreach_statement = 200,
  sym_foreach_parameters = 201,
  sym_foreach_parameter = 202,
  sym_reduction_list = 203,
  sym_reduction = 204,
  sym_reduction_operator = 205,
  sym_reduction_array = 206,
  sym_foreach_inner_statement = 207,
  sym_foreach_dimension_statement = 208,
  sym_forin_declaration_statement = 209,
  sym_forin_statement = 210,
  sym_forin_arguments = 211,
  sym_event_definition = 212,
  sym_event_parameters = 213,
  sym_event_parameter = 214,
  sym_boundary_definition = 215,
  sym_external_foreach_dimension = 216,
  sym_attribute = 217,
  sym_new_field = 218,
  sym_root = 219
};
#line 609 "/home/esteban/ProgramFile/basilisk/src/grid/config.h"

enum typedef_kind_t {
  sym_SCALAR = sym_root + 1,
  sym_VECTOR,
  sym_TENSOR,
  sym_COORD,
  sym__COORD,
  sym_VEC4,
  sym_IVEC
};
#line 21 "/home/esteban/ProgramFile/basilisk/src/common.h"


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;




int N = 16;


typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;


  scalar z;

} vector;

typedef struct {
  scalar * x;

  scalar * y;


  scalar * z;

} vectorl;

typedef struct {
  vector x;

  vector y;


  vector z;

} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
      omp_out.x += omp_in.x,
      omp_out.y += omp_in.y,
      omp_out.z += omp_in.z))
#line 95 "/home/esteban/ProgramFile/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  
    norm += sq(n->x);
    
#line 99
norm += sq(n->y);
    
#line 99
norm += sq(n->z);
  norm = sqrt(norm);
  
    n->x /= norm;
    
#line 102
n->y /= norm;
    
#line 102
n->z /= norm;
}

void origin (double x, double y, double z) {
  X0 = x; Y0 = y; Z0 = z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }
#line 122 "/home/esteban/ProgramFile/basilisk/src/common.h"
  enum { right, left, top, bottom, front, back };

int nboundary = 2*3;



#define _dirichlet(expr, ...) (2.*(expr) - val(_s,0,0,0))
#define _dirichlet_homogeneous(...) (- val(_s,0,0,0))
#define _dirichlet_face(expr,...) (expr)
#define _dirichlet_face_homogeneous(...) (0.)
#define _neumann(expr,...) (Delta*(expr) + val(_s,0,0,0))
#define _neumann_homogeneous(...) (val(_s,0,0,0))

double * _constant = NULL;
size_t datasize = 0;
typedef struct _Point Point;

#line 1 "grid/boundaries.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/esteban/ProgramFile/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 140 "/home/esteban/ProgramFile/basilisk/src/common.h"



typedef struct {
  int x;

  int y;


  int z;

} ivec;
typedef double (* BoundaryFunc) (Point, Point, scalar, bool *);
typedef struct {
  BoundaryFunc * boundary;
  BoundaryFunc * boundary_homogeneous;
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  ivec d;
  vector v;
  int face;
  bool nodump, freed;
  int block;
  scalar * depends;

  
#line 19 "/home/esteban/ProgramFile/basilisk/src/grid/stencils.h"
bool input, output;
  int width;
  int dirty;
  
#line 18 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid-common.h"
void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);
  
#line 28 "/home/esteban/ProgramFile/basilisk/src/vof.h"
scalar * tracers, c;
  bool inverse;
  
#line 21 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
scalar phi;
  
#line 460 "/home/esteban/ProgramFile/basilisk/src/heights.h"
vector height;
  
#line 22 "/home/esteban/ProgramFile/basilisk/src/tension.h"
double sigma;

#line 165 "/home/esteban/ProgramFile/basilisk/src/common.h"
} _Attributes;

static _Attributes * _attribute = NULL;






int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ ns++;}}
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s1=*_i;(&s1)->i>=0;s1=*++_i){
      if (s1.i == s.i)
 return true;}}
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      list = list_append (list, s);}}
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  {scalar*_i=(scalar*)( l2);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    l3 = list_append (l3, s);}}
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);}}
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){ nv++;}}
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  {vector*_i=(vector*)( list);if(_i)for(vector w=*_i;(&w)->x.i>=0;w=*++_i){ {
    bool id = true;
    
      if (w.x.i != v.x.i)
 id = false;
      
#line 266
if (w.y.i != v.y.i)
 id = false;
      
#line 266
if (w.z.i != v.z.i)
 id = false;
    if (id)
      return list;
  }}}
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    {vector*_i=(vector*)( l);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
      list = vectors_append (list, v);}}
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
     {
      if (!(s->i >= 0)) qassert ("/home/esteban/ProgramFile/basilisk/src/common.h", 289, "s->i >= 0");
      v.x = *s++;
    } 
#line 288
{
      if (!(s->i >= 0)) qassert ("/home/esteban/ProgramFile/basilisk/src/common.h", 289, "s->i >= 0");
      v.y = *s++;
    } 
#line 288
{
      if (!(s->i >= 0)) qassert ("/home/esteban/ProgramFile/basilisk/src/common.h", 289, "s->i >= 0");
      v.z = *s++;
    }
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  {tensor*_i=(tensor*)( list);if(_i)for(tensor t=*_i;(&t)->x.x.i>=0;t=*++_i){ nt++;}}
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
     {
      if (!(v->x.i >= 0)) qassert ("/home/esteban/ProgramFile/basilisk/src/common.h", 320, "v->x.i >= 0");
      t.x = *v++;
    } 
#line 319
{
      if (!(v->y.i >= 0)) qassert ("/home/esteban/ProgramFile/basilisk/src/common.h", 320, "v->x.i >= 0");
      t.y = *v++;
    } 
#line 319
{
      if (!(v->z.i >= 0)) qassert ("/home/esteban/ProgramFile/basilisk/src/common.h", 320, "v->x.i >= 0");
      t.z = *v++;
    }
    list = tensors_append (list, t);
  }
  return list;
}

static inline bool is_vertex_scalar (scalar s)
{
  
    if (_attribute[s.i].d.x != -1)
      return false;
    
#line 331
if (_attribute[s.i].d.y != -1)
      return false;
    
#line 331
if (_attribute[s.i].d.z != -1)
      return false;
  return true;
}

scalar * all = NULL;
scalar * baseblock = NULL;



scalar (* init_scalar) (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
vector (* init_face_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
void (* scalar_clone) (scalar, scalar);






static double mpi_time = 0.;


typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);

  t.tm = mpi_time;

  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



const vector zerof = {{_NVARMAX+4},{_NVARMAX+5},{_NVARMAX+6}};
const vector unityf = {{_NVARMAX+7},{_NVARMAX+8},{_NVARMAX+9}};
const scalar unity = {_NVARMAX+10};
const scalar zeroc = {_NVARMAX+11};



        vector fm = {{_NVARMAX+12},{_NVARMAX+13},{_NVARMAX+14}};
        scalar cm = {_NVARMAX+15};
#line 407 "/home/esteban/ProgramFile/basilisk/src/common.h"
void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = 1e30;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 do { double _tmp_ = m[irow][l]; m[irow][l] = m[icol][l]; m[icol][l] = _tmp_; } while(false);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 do { double _tmp_ = m[k][indxr[l]]; m[k][indxr[l]] = m[k][indxc[l]]; m[k][indxc[l]] = _tmp_; } while(false);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}



typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}



static char * display_defaults = NULL;

static void free_display_defaults() {
  pfree (display_defaults,__func__,__FILE__,__LINE__);
}

void display (const char * commands, bool overwrite)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (overwrite) {
    pfree (display_defaults,__func__,__FILE__,__LINE__);
    display_defaults = pmalloc (strlen(commands) + 2,__func__,__FILE__,__LINE__);
    strcpy (display_defaults, "@");
    strcat (display_defaults, commands);
  }
  else {
    if (!display_defaults)
      display_defaults = pstrdup ("@",__func__,__FILE__,__LINE__);
    display_defaults =
      prealloc (display_defaults,
        strlen(display_defaults) + strlen(commands) + 1,__func__,__FILE__,__LINE__);
    strcat (display_defaults, commands);
  }
}



typedef struct {
  double x;

  double y;


  double z;

} _coord;



typedef struct {
  float r, g, b, a;
} vec4;

#define attroffset(x) (offsetof(_Attributes,x))
#line 542 "/home/esteban/ProgramFile/basilisk/src/common.h"
#define BEGIN_FOREACH
#define END_FOREACH





#line 1 "grid/stencils.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/stencils.h"
#line 17 "/home/esteban/ProgramFile/basilisk/src/grid/stencils.h"










typedef struct _External External;

struct _External {
  char * name;
  void * pointer;
  int type;
  int nd;
  char reduct;
  char global;
  void * data;
  scalar s;
  External * externals, * next;
  int used;
};

typedef struct {
  const char * fname;
  int line;
  int first;
  int face;
  bool vertex;
  int parallel;
  scalar * listc;
  vectorl listf;
  scalar * dirty;
  void * data;
} ForeachData;


#define foreach_stencil(...) {\
  static int _first = 1.;\
  ForeachData _loop = {\
    .fname = __FILE__, .line = __LINE__, .first = _first\
  };\
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock;\
  s.i >= 0; i++, s = *i) {\
    _attribute[s.i].input = _attribute[s.i].output = false;\
    _attribute[s.i].width = 0;\
  }\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0}; NOT_UNUSED (point);\

#line 68


#define end_foreach_stencil()\
  check_stencil (&_loop);\
  boundary_stencil (&_loop);\
  _first = 0;\
}\

#line 75


#define foreach_vertex_stencil(...) foreach_stencil(__VA_ARGS__) _loop.vertex = true;
#define end_foreach_vertex_stencil() end_foreach_stencil()

#define foreach_face_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_face_stencil() end_foreach_stencil()

#define foreach_visible_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_visible_stencil() end_foreach_stencil()

#define foreach_level_stencil(...) {\
  if (0) {\
\
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
    Point point = {0}; NOT_UNUSED (point);\

#line 91

#define end_foreach_level_stencil() }}

#define foreach_coarse_level_stencil(...) foreach_level_stencil(__VA_ARGS__)
#define end_foreach_coarse_level_stencil() end_foreach_level_stencil()

#define foreach_level_or_leaf_stencil(...) foreach_level_stencil(__VA_ARGS__)
#define end_foreach_level_or_leaf_stencil() end_foreach_level_stencil()

#define foreach_point_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_point_stencil() end_foreach_stencil()

#define foreach_region_stencil(...) foreach_stencil(__VA_ARGS__)
#define end_foreach_region_stencil() end_foreach_stencil()

#define _stencil_is_face_x() { _loop.face |= (1 << 0);
#define end__stencil_is_face_x() }
#define _stencil_is_face_y() { _loop.face |= (1 << 1);
#define end__stencil_is_face_y() }
#define _stencil_is_face_z() { _loop.face |= (1 << 2);
#define end__stencil_is_face_z() }

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow);
void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line);

#define _stencil_val(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, false)\

#line 120

#define _stencil_val_o(a,_i,_j,_k)\
  stencil_val (point, a, _i, _j, _k, __FILE__, __LINE__, true)\

#line 123

#define _stencil_val_a(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, false, __FILE__, __LINE__)\

#line 126

#define _stencil_val_r(a,_i,_j,_k)\
  stencil_val_a (point, a, _i, _j, _k, true, __FILE__, __LINE__)\

#line 129


#define _stencil_fine(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_fine(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_fine_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_fine_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define _stencil_coarse(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
#define _stencil_coarse_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
#define _stencil_coarse_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

#define r_assign(x)
#define _assign(x)

#define _stencil_neighbor(i,j,k)
#define _stencil_child(i,j,k)
#define _stencil_aparent(i,j,k)
#define _stencil_aparent_a(i,j,k)
#define _stencil_aparent_r(i,j,k)

#define _stencil_neighborp(i,j,k) neighborp(i,j,k)

int _stencil_nop;
#define _stencil_val_higher_dimension (_stencil_nop = 1)
#define _stencil__val_constant(a,_i,_j,_k) (_stencil_nop = 1)
#define _stencil_val_diagonal(a,_i,_j,_k) (_stencil_nop = 1)

typedef void _stencil_undefined;

#define o_stencil -2







static inline bool scalar_is_dirty (scalar s)
{
  if (_attribute[s.i].dirty)
    return true;
  scalar * depends = _attribute[s.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      return true;}}
  return false;
}




static inline bool scalar_depends_from (scalar a, scalar b)
{
  scalar * depends = _attribute[a.i].depends;
  {scalar*_i=(scalar*)( depends);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (s.i == b.i)
      return true;}}
  return false;
}







void boundary_internal (scalar * list, const char * fname, int line);
void (* boundary_face) (vectorl);







void check_stencil (ForeachData * loop)
{
  loop->listf = (vectorl){NULL};




  {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    bool write = _attribute[s.i].output, read = _attribute[s.i].input;




    {





      if (read && scalar_is_dirty (s)) {





 if (_attribute[s.i].face) {
   if (_attribute[s.i].width > 0)
     loop->listc = list_append (loop->listc, s);
   else if (!write) {
     scalar sn = _attribute[s.i].v.x.i >= 0 ? _attribute[s.i].v.x : s;
     
       if (_attribute[s.i].v.x.i == s.i) {




  if (_attribute[sn.i].boundary[left] || _attribute[sn.i].boundary[right])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.x = list_append (loop->listf.x, s);
       }
       
#line 235
if (_attribute[s.i].v.y.i == s.i) {




  if (_attribute[sn.i].boundary[bottom] || _attribute[sn.i].boundary[top])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.y = list_append (loop->listf.y, s);
       }
       
#line 235
if (_attribute[s.i].v.z.i == s.i) {




  if (_attribute[sn.i].boundary[back] || _attribute[sn.i].boundary[front])
    loop->listc = list_append (loop->listc, s);
  else if (_attribute[s.i].dirty != 2)
    loop->listf.z = list_append (loop->listf.z, s);
       }
   }
 }





 else if (_attribute[s.i].width > 0)
   loop->listc = list_append (loop->listc, s);
      }





      if (write) {
 if (3 > 1 && !loop->vertex && loop->first) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 264
if (_attribute[s.i].d.y != -1)
       vertex = false;
     
#line 264
if (_attribute[s.i].d.z != -1)
       vertex = false;
   if (vertex)
     fprintf (ferr,
       "%s:%d: warning: vertex scalar '%s' should be assigned with"
       " a foreach_vertex() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 if (_attribute[s.i].face) {
   if (loop->face == 0 && loop->first)
     fprintf (ferr,
       "%s:%d: warning: face vector '%s' should be assigned with"
       " a foreach_face() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 else if (loop->face) {
   if (_attribute[s.i].v.x.i < 0) {
     int d = 1, i = 0;
      {
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.x.i = s.i;
  _attribute[s.i].boundary[left] = _attribute[s.i].boundary[right] = NULL;





       }
       d *= 2, i++;
     } 
#line 282
{
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.y.i = s.i;
  _attribute[s.i].boundary[bottom] = _attribute[s.i].boundary[top] = NULL;





       }
       d *= 2, i++;
     } 
#line 282
{
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.z.i = s.i;
  _attribute[s.i].boundary[back] = _attribute[s.i].boundary[front] = NULL;





       }
       d *= 2, i++;
     }
     if (!_attribute[s.i].face && loop->first)
       fprintf (ferr,
         "%s:%d: warning: scalar '%s' should be assigned with "
         "a foreach_face(x|y|z) loop\n",
         loop->fname, loop->line, _attribute[s.i].name);
   }
   else {
     char * name = NULL;
     if (_attribute[s.i].name) {
       name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
       char * s = name + strlen(name) - 1;
       while (s != name && *s != '.') s--;
       if (s != name) *s = '\0';
     }
     struct { int x, y, z; } input, output;
     vector v = _attribute[s.i].v;

     
       input.x = _attribute[v.x.i].input, output.x = _attribute[v.x.i].output;
       
#line 312
input.y = _attribute[v.y.i].input, output.y = _attribute[v.y.i].output;
       
#line 312
input.z = _attribute[v.z.i].input, output.z = _attribute[v.z.i].output;

     init_face_vector (v, name);


     
       _attribute[v.x.i].input = input.x, _attribute[v.x.i].output = output.x;
       
#line 318
_attribute[v.y.i].input = input.y, _attribute[v.y.i].output = output.y;
       
#line 318
_attribute[v.z.i].input = input.z, _attribute[v.z.i].output = output.z;





     pfree (name,__func__,__FILE__,__LINE__);
   }
 }
 else if (loop->vertex) {
   bool vertex = true;
   
     if (_attribute[s.i].d.x != -1)
       vertex = false;
     
#line 330
if (_attribute[s.i].d.y != -1)
       vertex = false;
     
#line 330
if (_attribute[s.i].d.z != -1)
       vertex = false;
   if (!vertex) {
     char * name = NULL;
     if (_attribute[s.i].name) name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     init_vertex_scalar (s, name);
     
       _attribute[s.i].v.x.i = -1;
       
#line 337
_attribute[s.i].v.y.i = -1;
       
#line 337
_attribute[s.i].v.z.i = -1;




     pfree (name,__func__,__FILE__,__LINE__);
   }
 }





 loop->dirty = list_append (loop->dirty, s);
 {scalar*_i=(scalar*)( baseblock);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
   if (scalar_depends_from (d, s))
     loop->dirty = list_append (loop->dirty, d);}}
      }
    }
  }}}
}




void boundary_stencil (ForeachData * loop)
{
  bool flux = false;
  
    if (loop->listf.x)
      flux = true;
    
#line 366
if (loop->listf.y)
      flux = true;
    
#line 366
if (loop->listf.z)
      flux = true;
  if (flux) {
#line 381 "/home/esteban/ProgramFile/basilisk/src/grid/stencils.h"
    boundary_face (loop->listf);
    
      pfree (loop->listf.x,__func__,__FILE__,__LINE__), loop->listf.x = NULL;
      
#line 383
pfree (loop->listf.y,__func__,__FILE__,__LINE__), loop->listf.y = NULL;
      
#line 383
pfree (loop->listf.z,__func__,__FILE__,__LINE__), loop->listf.z = NULL;
  }




  if (loop->listc) {






    boundary_internal (loop->listc, loop->fname, loop->line);
    pfree (loop->listc,__func__,__FILE__,__LINE__), loop->listc = NULL;
  }





  if (loop->dirty) {






    {scalar*_i=(scalar*)( loop->dirty);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = true;}}
    pfree (loop->dirty,__func__,__FILE__,__LINE__), loop->dirty = NULL;
  }
}
#line 550 "/home/esteban/ProgramFile/basilisk/src/common.h"





typedef struct {
  coord x, y, z;
} mat3;

OMP(omp declare reduction (+ : mat3 :
      omp_out.x.x += omp_in.x.x,
      omp_out.x.y += omp_in.x.y,
      omp_out.x.z += omp_in.x.z,
      omp_out.y.x += omp_in.y.x,
      omp_out.y.y += omp_in.y.y,
      omp_out.y.z += omp_in.y.z,
      omp_out.z.x += omp_in.z.x,
      omp_out.z.y += omp_in.z.y,
      omp_out.z.z += omp_in.z.z
      ))
#line 14 "main3D-cpp.c"
#line 1 "main3D.c"
#line 9 "main3D.c"
#line 1 "grid/multigrid3D.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid3D.h"

#line 1 "grid/multigrid.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"



typedef double real;
#line 24 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"
typedef struct {
  Grid g;
  char * d;
  size_t field_size;
} Multigrid;

struct _Point {
  int i;

  int j;


  int k;

  int level, n;
#ifdef foreach_block
  int l;
  #define _BLOCK_INDEX , point.l
#else
  #define _BLOCK_INDEX
#endif
};
static Point last_point;
#line 77 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"
#undef val
#define val(a,l,m,o) (((real *)((Multigrid *)grid)->d)[point.k + (o) +\
       ((1 << point.level) + 2*2)*\
       (point.j + (m) +\
        (point.i + (l))*((1 << point.level) + 2*2)) +\
       (((1 << 3*(point.level)) - 1)/7 + 2*2*((1 << 2*(point.level)) - 1) + 12*sq(2)*((1 << (point.level)) - 1) + 8*cube(2)*(point.level)) +\
       _index(a,0)*((Multigrid *)grid)->field_size])\

#line 84

#line 113 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"
#define allocated(a,l,m) (point.i+(a) >= 0 &&\
         point.i+(a) < (1 << point.level) + 2*2 &&\
         point.j+(l) >= 0 &&\
         point.j+(l) < (1 << point.level) + 2*2 &&\
         point.k+(m) >= 0 &&\
         point.k+(m) < (1 << point.level) + 2*2)\

#line 119


#define allocated_child(a,l,m) (level < depth() &&\
         point.i > 0 && point.i <= (1 << point.level) + 2 &&\
         point.j > 0 && point.j <= (1 << point.level) + 2 &&\
         point.k > 0 && point.k <= (1 << point.level) + 2)\

#line 125




#define depth() (grid->depth)
#line 173 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"
#define fine(a,l,m,o)\
(((real *)((Multigrid *)grid)->d)[2*point.k - 2 + (o) +\
   ((1 << point.level)*2 + 2*2)*\
   (2*point.j - 2 + (m) +\
    (2*point.i - 2 + (l))*((1 << point.level)*2 + 2*2)) +\
   (((1 << 3*(point.level + 1)) - 1)/7 + 2*2*((1 << 2*(point.level + 1)) - 1) + 12*sq(2)*((1 << (point.level + 1)) - 1) + 8*cube(2)*(point.level + 1)) +\
   _index(a,0)*((Multigrid *)grid)->field_size])\

#line 180

#define coarse(a,l,m,o)\
(((real *)((Multigrid *)grid)->d)[(point.k + 2)/2 + (o) +\
   ((1 << point.level)/2 + 2*2)*\
   ((point.j + 2)/2 + (m) +\
    ((point.i + 2)/2 + (l))*((1 << point.level)/2 + 2*2)) +\
   (((1 << 3*(point.level - 1)) - 1)/7 + 2*2*((1 << 2*(point.level - 1)) - 1) + 12*sq(2)*((1 << (point.level - 1)) - 1) + 8*cube(2)*(point.level - 1)) +\
   _index(a,0)*((Multigrid *)grid)->field_size])\

#line 188

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
  struct { int x, y, z; } child = {\
    2*((point.i + 2)%2) - 1,\
    2*((point.j + 2)%2) - 1,\
    2*((point.k + 2)%2) - 1\
  }; NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2;\
  parent.j = (point.j + 2)/2;\
  parent.k = (point.k + 2)/2;\

#line 202



#define foreach_level(l)\
OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k < point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j < point.n + 2; point.j++)\
\
      for (point.k = 2; point.k < point.n + 2; point.k++)\
\
 {\
\
          POINT_VARIABLES\

#line 222

#define end_foreach_level()\
\
 }\
\
  }\
}\

#line 229


#define foreach()\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = depth(); point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k < point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j < point.n + 2; point.j++)\
\
      for (point.k = 2; point.k < point.n + 2; point.k++)\
\
 {\
\
          POINT_VARIABLES\

#line 248

#define end_foreach()\
\
 }\
\
  }\
}\

#line 255


#define is_active(cell) (true)
#define is_leaf(cell) (level == depth())
#define is_local(cell) (true)
#define leaf 2
#define refine_cell(...) do {\
  fprintf (stderr, "grid depths do not match. Aborting.\n");\
  if (!(0)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h", 263, "0");\
} while (0)\

#line 265

#define tree ((Multigrid *)grid)
#line 1 "grid/foreach_cell.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/esteban/ProgramFile/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 4: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 5: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 6: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 7: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
\
\
  Point root = {2,2,2,0};\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = {0};\
  for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)\
\
    for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)\
\
\
      for (root.k = 2*Period.z; root.k <= 2*(2 - Period.z); root.k++)\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 4:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 5:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 6:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 7:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
\
\
    Point root = {2,2,2,0};\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = {0};\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
      for (root.k = 0; root.k <= 2*2; root.k++)\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 268 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"

#define foreach_face_generic()\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = depth(); point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k <= point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j <= point.n + 2; point.j++)\
\
      for (point.k = 2; point.k <= point.n + 2; point.k++)\
\
        {\
\
   POINT_VARIABLES\

#line 286

#define end_foreach_face_generic()\
\
 }\
\
  }\
}\

#line 293


#define foreach_vertex()\
foreach_face_generic() {\
  x -= Delta/2.;\
\
  y -= Delta/2.;\
\
\
  z -= Delta/2.;\
\

#line 304

#define end_foreach_vertex() } end_foreach_face_generic()

#define is_coarse() (point.level < depth())
#line 359 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"
#define foreach_vertex_aux()\
foreach_vertex() {\
  struct { int x, y, z; } _a = {point.i, point.j, point.k};\

#line 362

#define end_foreach_vertex_aux() } end_foreach_vertex()






#define is_face_x() { int ig = -1; VARIABLES; if (point.j < point.n + 2 && point.k < point.n + 2) {
#define end_is_face_x() }}
#define is_face_y() { int jg = -1; VARIABLES; if (point.i < point.n + 2 && point.k < point.n + 2) {
#define end_is_face_y() }}
#define is_face_z() { int kg = -1; VARIABLES; if (point.i < point.n + 2 && point.j < point.n + 2) {
#define end_is_face_z() }}

#define foreach_child() {\
  int _i = 2*point.i - 2;\
  int _j = 2*point.j - 2;\
  int _k = 2*point.k - 2;\
  point.level++;\
  point.n *= 2;\
  for (int _l = 0; _l < 2; _l++)\
    for (int _m = 0; _m < 2; _m++)\
      for (int _n = 0; _n < 2; _n++) {\
 point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;\
 POINT_VARIABLES;\

#line 388

#define end_foreach_child()\
  }\
  point.i = (_i + 2)/2;\
  point.j = (_j + 2)/2;\
  point.k = (_k + 2)/2;\
  point.level--;\
  point.n /= 2;\
}\

#line 397

#define foreach_child_break() _l = _m = _n = 2


#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

#line 1 "grid/neighbors.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/neighbors.h"
#line 35 "/home/esteban/ProgramFile/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j, _k = point.k;\
  for (int _l = - _nn; _l <= _nn; _l++) {\
    point.i = _i + _l;\
    for (int _m = - _nn; _m <= _nn; _m++) {\
      point.j = _j + _m;\
      for (int _n = - _nn; _n <= _nn; _n++) {\
 point.k = _k + _n;\
 POINT_VARIABLES;\

#line 45

#define end_foreach_neighbor()\
      }\
    }\
  }\
  point.i = _i; point.j = _j; point.k = _k;\
}\

#line 52

#define foreach_neighbor_break() _l = _m = _n = _nn + 1
#line 407 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s))
      for (int b = 0; b < _attribute[s.i].block; b++) {
 real * data = (real *) ((Multigrid *)grid)->d;
 data += (s.i + b)*((Multigrid *)grid)->field_size;
 for (int i = 0; i < ((Multigrid *)grid)->field_size; i++, data++)
   *data = val;
      }}}
}
#line 492 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"
#define foreach_boundary_dir(l,d)\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l < 0 ? depth() : l;\
  point.n = 1 << point.level;\
  int * _i = &point.j, * _j = &point.k;\
  if (d == left) {\
    point.i = 2;\
    ig = -1;\
  }\
  else if (d == right) {\
    point.i = point.n + 2 - 1;\
    ig = 1;\
  }\
  else if (d == bottom) {\
    point.j = 2;\
    _i = &point.i;\
    jg = -1;\
  }\
  else if (d == top) {\
    point.j = point.n + 2 - 1;\
    _i = &point.i;\
    jg = 1;\
  }\
  else if (d == back) {\
    point.k = 2;\
    _i = &point.i; _j = &point.j;\
    kg = -1;\
  }\
  else if (d == front) {\
    point.k = point.n + 2 - 1;\
    _i = &point.i; _j = &point.j;\
    kg = 1;\
  }\
  int _l;\
  OMP(omp for schedule(static))\
  for (_l = 0; _l < point.n + 2*2; _l++) {\
    *_i = _l;\
    for (int _m = 0; _m < point.n + 2*2; _m++) {\
      *_j = _m;\
      POINT_VARIABLES\

#line 534

#define end_foreach_boundary_dir()\
    }\
  }\
}\

#line 539


#define neighbor(o,p,q)\
  ((Point){point.i+o, point.j+p, point.k+q, point.level, point.n _BLOCK_INDEX})\

#line 543

#define is_boundary(point) (point.i < 2 || point.i >= point.n + 2 ||\
    point.j < 2 || point.j >= point.n + 2 ||\
    point.k < 2 || point.k >= point.n + 2)\

#line 547




#define foreach_boundary(b)\
  if (default_scalar_bc[b] != periodic_bc)\
    foreach_boundary_dir (depth(), b)\
      if (!is_boundary(point)) {\

#line 555

#define end_foreach_boundary() } end_foreach_boundary_dir()

#define neighborp(k,l,o) neighbor(k,l,o)

static double periodic_bc (Point point, Point neighbor, scalar s, bool * data);

static void box_boundary_level (const Boundary * b, scalar * scalars, int l)
{
  extern double (* default_scalar_bc[]) (Point, Point, scalar, bool *);
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  for (int d = 0; d < 2*3; d++)
    if (default_scalar_bc[d] == periodic_bc)
      {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (!is_constant(s) && _attribute[s.i].block > 0) {
   if (is_vertex_scalar (s))
     _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
   else if (_attribute[s.i].face) {
     vector v = _attribute[s.i].v;
     _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
   }
 }}}
  for (int bghost = 1; bghost <= 2; bghost++)
    for (int d = 0; d < 2*3; d++) {

      scalar * list = NULL, * listb = NULL;
      {scalar*_i=(scalar*)( scalars);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (!is_constant(s) && _attribute[s.i].block > 0) {
   scalar sb = s;

   if (_attribute[s.i].v.x.i >= 0) {

     int j = 0;
     while ((&_attribute[s.i].v.x)[j].i != s.i) j++;
     sb = (&_attribute[s.i].v.x)[(j - d/2 + 3) % 3];
   }

   if (_attribute[sb.i].boundary[d] && _attribute[sb.i].boundary[d] != periodic_bc) {
     list = list_append (list, s);
     listb = list_append (listb, sb);
   }
 }}}

      if (list) {
  BEGIN_FOREACH{foreach_boundary_dir (l, d) {
   scalar s, sb;
   {scalar*_i0=listb;scalar*_i1= list;if(_i0)for(sb=*_i0,s=*_i1;_i0->i>= 0;sb=*++_i0,s=*++_i1){ {
     if ((_attribute[s.i].face && sb.i == _attribute[s.i].v.x.i) || is_vertex_scalar (s)) {

       if (bghost == 1)
 
    val(s,(ig + 1)/2,(jg + 1)/2,(kg + 1)/2) =
    _attribute[sb.i].boundary[d] (point, neighborp(ig,jg,kg), s, NULL);
     }
     else

      
  val(s,bghost*ig,bghost*jg,bghost*kg) =
  _attribute[sb.i].boundary[d] (neighborp((1 - bghost)*ig,
       (1 - bghost)*jg,
       (1 - bghost)*kg),
    neighborp(bghost*ig,bghost*jg,bghost*kg),
    s, NULL);
   }}}
 }end_foreach_boundary_dir();}END_FOREACH 
 pfree (list,__func__,__FILE__,__LINE__);
 pfree (listb,__func__,__FILE__,__LINE__);
      }
    }
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#line 755 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"
void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Multigrid * m = ((Multigrid *)grid);
  pfree (m->d,__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
  grid = NULL;
}

int log_base2 (int n) {
  int m = n, r = 0;
  while (m > 1)
    m /= 2, r++;
  return (1 << r) < n ? r + 1 : r;
}

void init_grid (int n)
{
  free_grid();
  Multigrid * m = ((Multigrid *) pmalloc ((1)*sizeof(Multigrid),__func__,__FILE__,__LINE__));
  grid = (Grid *) m;
  grid->depth = grid->maxdepth = log_base2(n);
  m->field_size = (((1 << 3*(depth() + 1)) - 1)/7 + 2*2*((1 << 2*(depth() + 1)) - 1) + 12*sq(2)*((1 << (depth() + 1)) - 1) + 8*cube(2)*(depth() + 1));
  N = 1 << depth();

  grid->n = grid->tn = 1 << 3*depth();

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  add_boundary (b);

  Boundary * mpi_boundary_new();
  mpi_boundary_new();
#line 799 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"
  m->d = (char *) pmalloc(m->field_size*datasize,__func__,__FILE__,__LINE__);
  reset (all, 0.);
}

void realloc_scalar (int size)
{
  Multigrid * p = ((Multigrid *)grid);
  datasize += size;
  p->d = (char *) prealloc (p->d, (p->field_size*datasize)*sizeof(char),__func__,__FILE__,__LINE__);
}


int mpi_dims[3], mpi_coords[3];
#line 824 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"
Point locate (double xp, double yp, double zp)
{
  Point point = {0};int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  point.level = -1, point.n = 1 << depth();

  point.i = (xp - X0)/L0*point.n*mpi_dims[0] + 2 - mpi_coords[0]*point.n;
  if (point.i < 2 || point.i >= point.n + 2)
    return point;

  point.j = (yp - Y0)/L0*point.n*mpi_dims[0] + 2 - mpi_coords[1]*point.n;
  if (point.j < 2 || point.j >= point.n + 2)
    return point;


  point.k = (zp - Z0)/L0*point.n*mpi_dims[0] + 2 - mpi_coords[2]*point.n;
  if (point.k < 2 || point.k >= point.n + 2)
    return point;
#line 857 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"
  point.level = depth();
  return point;
}


#line 1 "grid/multigrid-common.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/events.h"
typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
static void _init_solver (void);





static int END_EVENT = 1234567890;
static double TEND_EVENT = 1234567890;
static double TEPS = 1e-9;

static void event_error (Event * ev, const char * s)
{
  fprintf (ferr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = -1; ev->t = - TEND_EVENT;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = - TEND_EVENT;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == - TEND_EVENT) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = -1; ev->t = - TEND_EVENT;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == END_EVENT || ev->t == TEND_EVENT) {
 ev->i = END_EVENT; ev->t = - TEND_EVENT;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != - TEND_EVENT)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->i = -1; ev->t = - TEND_EVENT;
  return event_done;
}

void event_register (Event event) {
  if (!(Events)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/events.h", 114, "Events");
  if (!(!event.last)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/events.h", 115, "!event.last");
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      if (!(parent < 0)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/events.h", 119, "parent < 0");
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 162 "/home/esteban/ProgramFile/basilisk/src/grid/events.h"
static bool overload_event() { return true; }

static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (!overload_event() || iter == ev->i || fabs (t - ev->t) <= TEPS*t) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == END_EVENT && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = END_EVENT; tnext = 1e30;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != END_EVENT &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != END_EVENT &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if (overload_event() && (!cond || cond1) && (tnext != 1e30 || inext != END_EVENT)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != 1e30 && tnext > t) {
    if (!(dt > 0.)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/events.h", 277, "dt > 0.");
    unsigned int n = (tnext - t)/dt;
    if (!(n < INT_MAX)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/events.h", 279, "n < INT_MAX");
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt*(1. + TEPS))
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}

void init_solver()
{
  Events = pmalloc (sizeof (Event),__func__,__FILE__,__LINE__);
  Events[0].last = 1;
  _attribute = pcalloc (datasize/sizeof(real), sizeof (_Attributes),__func__,__FILE__,__LINE__);
  int n = datasize/sizeof(real);
  all = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  baseblock = (scalar *) pmalloc (sizeof (scalar)*(n + 1),__func__,__FILE__,__LINE__);
  for (int i = 0; i < n; i++)
    baseblock[i].i = all[i].i = i;
  baseblock[n].i = all[n].i = -1;




  mpi_init();





}
#line 2 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])
#define diagonalize(a)
#define val_diagonal(a,k,l,m) ((k) == 0 && (l) == 0 && (m) == 0)

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level)/mpi_dims[0]);\
  double Delta_x = Delta;\
\
  double Delta_y = Delta;\
\
\
  double Delta_z = Delta;\
\
\
  double x = (ig/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
  double z = (kg/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)) + 0.5)*Delta + Z0;\
\
\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  NOT_UNUSED(Delta_x);\
\
  NOT_UNUSED(Delta_y);\
\
\
  NOT_UNUSED(Delta_z);\
\
\
  ;\

#line 44


#line 1 "grid/fpe.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static int gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', ferr);
    fflush (ferr);
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  return system (command);
}

static void caught_abort (int sig)
{
  fprintf (ferr, "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (ferr, "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (ferr, "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 47 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()

#define foreach_point(...)\
{\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  coord _p = { __VA_ARGS__ };\
  Point point = locate (_p.x, _p.y, _p.z);\
  if (point.level >= 0) {\
    POINT_VARIABLES\

#line 57

#define end_foreach_point() }}

#define foreach_region(p, box, n)\
 OMP_PARALLEL() { NOT_UNUSED (p);\
    coord p = {0, 0, box[0].z};\
    OMP(omp for schedule(static))\
      for (int _i = 0; _i < (int) n.x; _i++) {\
 p.x = box[0].x + (box[1].x - box[0].x)/n.x*(_i + 0.5);\
 for (int _j = 0; _j < (int) n.y; _j++) {\
   p.y = box[0].y + (box[1].y - box[0].y)/n.y*(_j + 0.5);\
   Point point = locate (p.x, p.y, p.z);\
   if (point.level >= 0) {\
     int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
     POINT_VARIABLES\

#line 72

#define end_foreach_region() }}}}
#line 133 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
static void init_block_scalar (scalar sb, const char * name, const char * ext,
          int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    strcat (strcpy (bname, name), ext);
    _attribute[sb.i].block = block;
    baseblock = list_append (baseblock, sb);
  }
  else {

    sprintf (bname, "%s%d%s", name, n, ext);
    _attribute[sb.i].block = - n;
  }
  _attribute[sb.i].name = pstrdup (bname,__func__,__FILE__,__LINE__);
  all = list_append (all, sb);
}

#define interpreter_set_int(...)
#define interpreter_reset_scalar(...)

scalar alloc_block_scalar (const char * name, const char * ext, int block)
{
  interpreter_set_int (&block);
  int nvar = datasize/sizeof(real);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && _attribute[sb.i].freed)
      n++, sb.i++;
    if (n >= block) {
      memset (&_attribute[s.i], 0, block*sizeof (_Attributes));
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++) {
 init_block_scalar (sb, name, ext, n, block);
 interpreter_reset_scalar (sb);
      }
      trash (((scalar []){s, {-1}}));
      return s;
    }
    s.i = sb.i + 1;
  }


  s = (scalar){nvar};
  if (!(nvar + block <= _NVARMAX)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h", 179, "nvar + block <= _NVARMAX");

  if (_attribute == NULL)
    _attribute = (_Attributes *) pcalloc (nvar + block + 1, sizeof (_Attributes),__func__,__FILE__,__LINE__);
  else
    _attribute = (_Attributes *)
      prealloc (_attribute, (nvar + block + 1)*sizeof (_Attributes),__func__,__FILE__,__LINE__);
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }

  realloc_scalar (block*sizeof(real));
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  scalar s = alloc_block_scalar (name, ext, block), sb;
  int n = 0;
  for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
    init_scalar (sb, NULL);
  return s;
}

scalar new_scalar (const char * name)
{
  return init_scalar (alloc_block_scalar (name, "", 1), NULL);
}

scalar new_vertex_scalar (const char * name)
{
  return init_vertex_scalar (alloc_block_scalar (name, "", 1), NULL);
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  
    v.x = alloc_block_scalar (name, ext.x, block);
    
#line 221
v.y = alloc_block_scalar (name, ext.y, block);
    
#line 221
v.z = alloc_block_scalar (name, ext.z, block);
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_face_vector (v, NULL);
  return v;
}

vector new_block_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 245
vb.y.i = v.y.i + i;
      
#line 245
vb.z.i = v.z.i + i;
    init_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 248
_attribute[vb.y.i].block = - i;
      
#line 248
_attribute[vb.z.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 251
_attribute[v.y.i].block = block;
    
#line 251
_attribute[v.z.i].block = block;
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    
      vb.x.i = v.x.i + i;
      
#line 261
vb.y.i = v.y.i + i;
      
#line 261
vb.z.i = v.z.i + i;
    init_face_vector (vb, NULL);
    
      _attribute[vb.x.i].block = - i;
      
#line 264
_attribute[vb.y.i].block = - i;
      
#line 264
_attribute[vb.z.i].block = - i;
  }
  
    _attribute[v.x.i].block = block;
    
#line 267
_attribute[v.y.i].block = block;
    
#line 267
_attribute[v.z.i].block = block;
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  tensor t;
   {
    strcat (strcpy (cname, name), ext.x);
    t.x = alloc_block_vector (cname, 1);
  } 
#line 276
{
    strcat (strcpy (cname, name), ext.y);
    t.y = alloc_block_vector (cname, 1);
  } 
#line 276
{
    strcat (strcpy (cname, name), ext.z);
    t.z = alloc_block_vector (cname, 1);
  }
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  struct { char * x, * y, * z; } ext = {".x.x", ".y.y", ".z.z"};
  tensor t;
  
    t.x.x = alloc_block_scalar (name, ext.x, 1);
    
#line 289
t.y.y = alloc_block_scalar (name, ext.y, 1);
    
#line 289
t.z.z = alloc_block_scalar (name, ext.z, 1);

    t.x.y = alloc_block_scalar (name, ".x.y", 1);
    t.y.x = t.x.y;


    t.x.z = alloc_block_scalar (name, ".x.z", 1);
    t.z.x = t.x.z;
    t.y.z = alloc_block_scalar (name, ".y.z", 1);
    t.z.y = t.y.z;




  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    _constant = (double *) prealloc (_constant, (s.i - _NVARMAX + 1)*sizeof(double),__func__,__FILE__,__LINE__);
    for (int i = nconst; i < s.i - _NVARMAX; i++)
      _constant[i] = 0.;
    nconst = s.i - _NVARMAX + 1;
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  
    init_const_scalar (v.x, name, *val++);
    
#line 330
init_const_scalar (v.y, name, *val++);
    
#line 330
init_const_scalar (v.z, name, *val++);
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  
    v.x.i = _NVARMAX + i++;
    
#line 337
v.y.i = _NVARMAX + i++;
    
#line 337
v.z.i = _NVARMAX + i++;
  init_const_vector (v, name, val);
  return v;
}

static void cartesian_scalar_clone (scalar clone, scalar src)
{
  char * cname = _attribute[clone.i].name;
  BoundaryFunc * boundary = _attribute[clone.i].boundary;
  BoundaryFunc * boundary_homogeneous = _attribute[clone.i].boundary_homogeneous;
  if (!(_attribute[src.i].block > 0 && _attribute[clone.i].block == _attribute[src.i].block)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h", 347, "src.block > 0 && clone.block == src.block");
  pfree (_attribute[clone.i].depends,__func__,__FILE__,__LINE__);
  _attribute[clone.i] = _attribute[src.i];
  _attribute[clone.i].name = cname;
  _attribute[clone.i].boundary = boundary;
  _attribute[clone.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[clone.i].boundary[i] = _attribute[src.i].boundary[i];
    _attribute[clone.i].boundary_homogeneous[i] = _attribute[src.i].boundary_homogeneous[i];
  }
  _attribute[clone.i].depends = list_copy (_attribute[src.i].depends);
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(real), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  {scalar*_i=(scalar*)( l);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    scalar c = _attribute[s.i].block > 1 ? new_block_scalar("c", "", _attribute[s.i].block) : new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }}}
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    {
      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
      
#line 374
if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];
      
#line 374
if (_attribute[s.i].v.z.i >= 0 && map[_attribute[s.i].v.z.i] >= 0)
 _attribute[s.i].v.z.i = map[_attribute[s.i].v.z.i];}}}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    for (int i = 0; i < _attribute[f.i].block; i++) {
      scalar fb = {f.i + i};
      if (_attribute[f.i].delete)
 _attribute[f.i].delete (fb);
      pfree (_attribute[fb.i].name,__func__,__FILE__,__LINE__); _attribute[fb.i].name = NULL;
      pfree (_attribute[fb.i].boundary,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary = NULL;
      pfree (_attribute[fb.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary_homogeneous = NULL;
      pfree (_attribute[fb.i].depends,__func__,__FILE__,__LINE__); _attribute[fb.i].depends = NULL;
      _attribute[fb.i].freed = true;
    }
  }}}

  if (list == all) {
    all[0].i = -1;
    baseblock[0].i = -1;
    return;
  }

  trash (list);
  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    if (_attribute[f.i].block > 0) {
      scalar * s;
      for (s = all; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[_attribute[f.i].block].i >= 0; s++)
   s[0] = s[_attribute[f.i].block];
 s->i = -1;
      }
      for (s = baseblock; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[1].i >= 0; s++)
   s[0] = s[1];
 s->i = -1;
      }
    }
  }}}
}

void free_solver()
{
  if (!(_val_higher_dimension == 0.)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h", 425, "_val_higher_dimension == 0.");

  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  pfree (baseblock,__func__,__FILE__,__LINE__); baseblock = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree (_attribute,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;




  free_grid();
  qpclose_all();
#if TRACE
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_face) (vectorl);




void boundary_flux (vector * list) __attribute__ ((deprecated));

void boundary_flux (vector * list)
{
  vectorl list1 = {NULL};
  {vector*_i=(vector*)( list);if(_i)for(vector v=*_i;(&v)->x.i>=0;v=*++_i){
    {
      list1.x = list_append (list1.x, v.x);
      
#line 481
list1.y = list_append (list1.y, v.y);
      
#line 481
list1.z = list_append (list1.z, v.z);}}}
  boundary_face (list1);
  
    pfree (list1.x,__func__,__FILE__,__LINE__);
    
#line 484
pfree (list1.y,__func__,__FILE__,__LINE__);
    
#line 484
pfree (list1.z,__func__,__FILE__,__LINE__);
}

static scalar * list_add_depends (scalar * list, scalar s)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
    if (t.i == s.i)
      return list;}}
  scalar * list1 = list;
  {scalar*_i=(scalar*)( _attribute[s.i].depends);if(_i)for(scalar d=*_i;(&d)->i>=0;d=*++_i){
    if (_attribute[d.i].dirty)
      list1 = list_add_depends (list1, d);}}
  return list_append (list1, s);
}

     
void boundary_internal (scalar * list, const char * fname, int line)
{tracing("boundary_internal","/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h",500);
  if (list == NULL)
    {end_tracing("boundary_internal","/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h",503);return;}
  scalar * listc = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (scalar_is_dirty (s)) {
 if (_attribute[s.i].face && _attribute[s.i].dirty != 2)
   {
     if (_attribute[s.i].v.x.i == s.i)
       listf.x = list_add (listf.x, s), flux = true;
     
#line 512
if (_attribute[s.i].v.y.i == s.i)
       listf.y = list_add (listf.y, s), flux = true;
     
#line 512
if (_attribute[s.i].v.z.i == s.i)
       listf.z = list_add (listf.z, s), flux = true;}
 if (!is_constant(cm) && _attribute[cm.i].dirty)
   listc = list_add_depends (listc, cm);
 if (_attribute[s.i].face != 2)
   listc = list_add_depends (listc, s);
      }




    }}}
  if (flux) {
    boundary_face (listf);
    
      pfree (listf.x,__func__,__FILE__,__LINE__);
      
#line 527
pfree (listf.y,__func__,__FILE__,__LINE__);
      
#line 527
pfree (listf.z,__func__,__FILE__,__LINE__);
  }
  if (listc) {






    boundary_level (listc, -1);
    {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = false;}}
    pfree (listc,__func__,__FILE__,__LINE__);
  }
end_tracing("boundary_internal","/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h",541);}

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_face (vectorl list)
{
  
    {scalar*_i=(scalar*)( list.x);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 551
{scalar*_i=(scalar*)( list.y);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
    
#line 551
{scalar*_i=(scalar*)( list.z);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].dirty = 2;}}
}

static double symmetry (Point point, Point neighbor, scalar s, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return -val(s,0,0,0);
}

BoundaryFunc default_scalar_bc[] = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  int block = _attribute[s.i].block;
  BoundaryFunc * boundary = _attribute[s.i].boundary;
  BoundaryFunc * boundary_homogeneous = _attribute[s.i].boundary_homogeneous;
  _attribute[s.i].name = pname;
  if (block < 0)
    _attribute[s.i].block = block;
  else
    _attribute[s.i].block = block > 0 ? block : 1;

  _attribute[s.i].boundary = boundary ? boundary : (BoundaryFunc *) pmalloc (nboundary*sizeof (BoundaryFunc),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = boundary_homogeneous ? boundary_homogeneous :
    (BoundaryFunc *) pmalloc (nboundary*sizeof (BoundaryFunc),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*3 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
   {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  } 
#line 595
{
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  } 
#line 595
{
    _attribute[s.i].d.z = 0;
    _attribute[s.i].v.z.i = -1;
  }
  _attribute[s.i].face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  cartesian_init_scalar (s, name);
  
    _attribute[s.i].d.x = -1;
    
#line 607
_attribute[s.i].d.y = -1;
    
#line 607
_attribute[s.i].d.z = -1;
  for (int d = 0; d < nboundary; d++)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
  return s;
}

BoundaryFunc default_vector_bc[] = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.x);
      cartesian_init_scalar (v.x, cname);
    }
    else
      cartesian_init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  } 
#line 622
{
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.y);
      cartesian_init_scalar (v.y, cname);
    }
    else
      cartesian_init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  } 
#line 622
{
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.z);
      cartesian_init_scalar (v.z, cname);
    }
    else
      cartesian_init_scalar (v.z, NULL);
    _attribute[v.z.i].v = v;
  }

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*3 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
   {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  } 
#line 642
{
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  } 
#line 642
{
    _attribute[v.z.i].d.z = -1;
    _attribute[v.z.i].face = true;
  }
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
   {
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.x);
      cartesian_init_vector (t.x, cname);
    }
    else
      cartesian_init_vector (t.x, NULL);
  } 
#line 654
{
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.y);
      cartesian_init_vector (t.y, cname);
    }
    else
      cartesian_init_vector (t.y, NULL);
  } 
#line 654
{
    if (name) {
      char cname[strlen(name) + 3];
      strcat (strcpy (cname, name), ext.z);
      cartesian_init_vector (t.z, cname);
    }
    else
      cartesian_init_vector (t.z, NULL);
  }
#line 678 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
    if (!(false)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h", 678, "false");

  return t;
}

void output_cells (FILE * fp, coord c, double size)
{
   BEGIN_FOREACH{foreach() {
    bool inside = true;
    coord o = {x,y,z};
    
      if (inside && size > 0. &&
   (o.x > c.x + size || o.x < c.x - size))
 inside = false;
      
#line 689
if (inside && size > 0. &&
   (o.y > c.y + size || o.y < c.y - size))
 inside = false;
      
#line 689
if (inside && size > 0. &&
   (o.z > c.z + size || o.z < c.z - size))
 inside = false;
    if (inside) {
      Delta /= 2.;
#line 704 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
      for (int i = -1; i <= 1; i += 2) {
 fprintf (fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n\n",
   x - Delta, y - Delta, z + i*Delta,
   x - Delta, y + Delta, z + i*Delta,
   x + Delta, y + Delta, z + i*Delta,
   x + Delta, y - Delta, z + i*Delta,
   x - Delta, y - Delta, z + i*Delta);
 for (int j = -1; j <= 1; j += 2)
   fprintf (fp, "%g %g %g\n%g %g %g\n\n",
     x + i*Delta, y + j*Delta, z - Delta,
     x + i*Delta, y + j*Delta, z + Delta);
      }

    }
  }end_foreach();}END_FOREACH 
  fflush (fp);
}
#line 729 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"







    "  splot '%s' w l lc 0, "
    "'%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 1"
           " title columnhead(4+4*v)",

    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp, (coord){x,y,z}, 4.*Delta);
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){





    fprintf (fp, "x y z %s ", _attribute[v.i].name);}}

  fputc ('\n', fp);
#line 810 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++)
 for (int m = -2; m <= 2; m++) {
   {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
     fprintf (fp, "%g %g %g ",
       x + k*Delta + _attribute[v.i].d.x*Delta/2.,
       y + l*Delta + _attribute[v.i].d.y*Delta/2.,
       z + m*Delta + _attribute[v.i].d.z*Delta/2.);
     if (allocated(k,l,m))
       fprintf (fp, "%g ", val(v,k,l,m));
     else
       fputs ("n/a ", fp);
   }}}
   fputc ('\n', fp);
 }

  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }}}
  fclose (fp);

  fprintf (ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vertex_scalar = cartesian_init_vertex_scalar;
  init_vector = cartesian_init_vector;
  init_face_vector = cartesian_init_face_vector;
  init_tensor = cartesian_init_tensor;
  boundary_level = cartesian_boundary_level;
  boundary_face = cartesian_boundary_face;
  scalar_clone = cartesian_scalar_clone;
  debug = cartesian_debug;
}

tensor init_symmetric_tensor (tensor t, const char * name)
{
  return init_tensor (t, name);
}

static double interpolate_linear (Point point, scalar v,
      double xp, double yp, double zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
#line 885 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
  x = (xp - x)/Delta - _attribute[v.i].d.x/2.;
  y = (yp - y)/Delta - _attribute[v.i].d.y/2.;
  z = (zp - z)/Delta - _attribute[v.i].d.z/2.;
  int i = sign(x), j = sign(y), k = sign(z);
  x = fabs(x); y = fabs(y); z = fabs(z);

  return (((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
    (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y)*(1. - z) +
   ((val(v,0,0,k)*(1. - x) + val(v,i,0,k)*x)*(1. - y) +
    (val(v,0,j,k)*(1. - x) + val(v,i,j,k)*x)*y)*z);

}


#line 867
static void _stencil_interpolate_linear (Point point, scalar v,
_stencil_undefined * xp,_stencil_undefined * yp,_stencil_undefined * zp)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;         
#line 885 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
        
        
        
  
          

_stencil_val(v,0,0,0);_stencil_val(v, o_stencil,0,0);
_stencil_val(v,0,o_stencil,0); _stencil_val(v,o_stencil,o_stencil,0);
_stencil_val(v,0,0,o_stencil); _stencil_val(v,o_stencil,0,o_stencil);
_stencil_val(v,0,o_stencil,o_stencil); _stencil_val(v,o_stencil,o_stencil,o_stencil);

  
#line 891
return           
              
    
    ;

}

     
double interpolate (scalar v, double xp, double yp, double zp,
      bool linear)
{tracing("interpolate","/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h",899);
  double val = 1e30;
  foreach_point_stencil (1,DEPAREN({xp, yp, zp}),)
    { _stencil_interpolate_linear (point, v, NULL, NULL, NULL); _stencil_val(v,0,0,0);    }end_foreach_point_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (min:val)){
#line 903
foreach_point (xp, yp, zp)
    val = linear ? interpolate_linear (point, v, xp, yp, zp) : val(v,0,0,0);end_foreach_point();mpi_all_reduce_array(&val,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
  
#line 905
{end_tracing("interpolate","/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h",905);return val;}
end_tracing("interpolate","/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h",906);}

     
void interpolate_array (scalar * list, coord * a, int n, double * v,
   bool linear)
{tracing("interpolate_array","/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h",909);
  int len = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    len++;}}
  for (int i = 0; i < n; i++) {
    double * w = v;
#line 926 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      *(w++) = 1e30;}}
    foreach_point_stencil (1,DEPAREN({a[i].x, a[i].y, a[i].z}),) {   
      
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 { _stencil_val(s,0,0,0); _stencil_interpolate_linear (point, s, NULL, NULL, NULL);    }}}
    }end_foreach_point_stencil()
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(min:v[:len])){
#line 928
foreach_point (a[i].x, a[i].y, a[i].z) {
      int j = 0;
      {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 v[j++] = !linear ? val(s,0,0,0) : interpolate_linear (point, s, a[i].x, a[i].y, a[i].z);}}
    }end_foreach_point();mpi_all_reduce_array(v,double,MPI_MIN,len);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 

    
#line 934
v = w;
  }
end_tracing("interpolate_array","/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h",936);}



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].boundary = (BoundaryFunc *) prealloc (_attribute[s.i].boundary, nboundary*sizeof (BoundaryFunc),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (BoundaryFunc *)
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (BoundaryFunc),__func__,__FILE__,__LINE__);
  }}}
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      
 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
 
#line 956
_attribute[v.z.i].boundary[b] = _attribute[v.z.i].boundary_homogeneous[b] = symmetry;
 
#line 956
_attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }}}
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return val(s,0,0,0);
}

static void periodic_boundary (int d)
{

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (is_vertex_scalar (s))
      _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
    else
      _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;}}

  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }}}

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{





    if (!(dir <= back)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h", 997, "dir <= back");


  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}


double getvalue (Point point, scalar s, int i, int j, int k)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return val(s,i,j,k);
}

void default_stencil (Point p, scalar * list)
{
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].input = true, _attribute[s.i].width = 2;}}
}




static void write_stencil_index (int * index)
{
  fprintf (qstderr(), "[%d", index[0]);
  for (int d = 1; d < 3; d++)
    fprintf (qstderr(), ",%d", index[d]);
  fputs ("]", qstderr());
}

void stencil_val (Point p, scalar s, int i, int j, int k,
    const char * file, int line, bool overflow)
{
  if (is_constant(s) || s.i < 0)
    return;
  int index[] = {i, j, k};
  for (int d = 0; d < 3; d++)
    index[d] += (&p.i)[d];
  bool central = true;
  for (int d = 0; d < 3; d++) {
    if (!overflow && (index[d] > 2 || index[d] < - 2)) {
      fprintf (qstderr(), "%s:%d: error: stencil overflow: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
    if (index[d] != 0)
      central = false;
  }
  if (central) {
    if (!_attribute[s.i].output)
      _attribute[s.i].input = true;
  }
  else {
    _attribute[s.i].input = true;
    int d = 0;
     {
      if ((!_attribute[s.i].face || _attribute[s.i].v.x.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    } 
#line 1057
{
      if ((!_attribute[s.i].face || _attribute[s.i].v.y.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    } 
#line 1057
{
      if ((!_attribute[s.i].face || _attribute[s.i].v.z.i != s.i) && abs(index[d]) > _attribute[s.i].width)
 _attribute[s.i].width = abs(index[d]);
      d++;
    }
  }
}

void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
      const char * file, int line)
{
  if (is_constant(s) || s.i < 0)
    abort();
  int index[] = {i, j, k};
  for (int d = 0; d < 3; d++)
    index[d] += (&p.i)[d];
  for (int d = 0; d < 3; d++)
    if (index[d] != 0) {
      fprintf (qstderr(), "%s:%d: error: illegal write: %s",
        file, line, _attribute[s.i].name);
      write_stencil_index (index);
      fprintf (qstderr(), "\n");
      fflush (qstderr());
      abort();
    }
  if (input && !_attribute[s.i].output)
    _attribute[s.i].input = true;
  _attribute[s.i].output = true;
}




#define dimensional(...)

#define show_dimension_internal(...)
#define display_value(...)
#define interpreter_verbosity(...)
#line 4 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double sum = 0.;
  {foreach_child()
    sum += val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 3);
}


#line 26
static void _stencil_restriction_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;   
  
  {foreach_child()
    { _stencil_val(s,0,0,0); }end_foreach_child()}
  _stencil_val_a(s,0,0,0);    
}

static inline void restriction_volume_average (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;

#line 35
if(!is_constant(cm)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double sum = 0.;
  {foreach_child()
    sum += val(cm,0,0,0)*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 3)/(val(cm,0,0,0) + 1e-30);
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 35
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double sum = 0.;
  {foreach_child()
    sum += _const_cm*val(s,0,0,0);end_foreach_child()}
  val(s,0,0,0) = sum/(1 << 3)/(_const_cm + 1e-30);
}}

#line 40
}

static inline void face_average (Point point, vector v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
   {







      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0) +
        fine(v.x,0,0,1) + fine(v.x,0,1,1))/4.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0) +
  fine(v.x,2,0,1) + fine(v.x,2,1,1))/4.;

  } 
#line 44
{







      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,0,0,1) +
        fine(v.y,1,0,0) + fine(v.y,1,0,1))/4.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,0,2,1) +
  fine(v.y,1,2,0) + fine(v.y,1,2,1))/4.;

  } 
#line 44
{







      val(v.z,0,0,0) = (fine(v.z,0,0,0) + fine(v.z,1,0,0) +
        fine(v.z,0,1,0) + fine(v.z,1,1,0))/4.;
      val(v.z,0,0,1) = (fine(v.z,0,0,2) + fine(v.z,1,0,2) +
  fine(v.z,0,1,2) + fine(v.z,1,1,2))/4.;

  }
}

static inline void restriction_face (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  face_average (point, _attribute[s.i].v);
}

static inline void restriction_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  for (int i = 0; i <= 1; i++) {
    val(s,i,0,0) = fine(s,2*i,0,0);

    val(s,i,1,0) = fine(s,2*i,2,0);


    for (int j = 0; j <= 1; j++)
      val(s,i,j,1) = fine(s,2*i,2*j,2);

  }
}

static inline void no_restriction (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;}

static inline void no_data (Point point, scalar s) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = 1e30;end_foreach_child()}
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar[]){s,{-1}}));
  for (int l = grid->maxdepth - 1; l >= 0; l--) {
    foreach_coarse_level_stencil (1,DEPAREN({l}),) {
      {foreach_child()
        { _stencil_val(s,0,0,0);_stencil_val_a(w,0,0,0); }end_foreach_child()}
default_stencil (
      
#line 93
point,((scalar[]){ s,{-1}}));
      {foreach_child() {  
         _stencil_val(s,0,0,0); 
_stencil_val(w,0,0,0);
        
#line 96
_stencil_val_a(s,0,0,0); 

        _stencil_val_r(w,0,0,0);  
      }end_foreach_child()}
    }end_foreach_coarse_level_stencil()
     BEGIN_FOREACH{
#line 90
foreach_coarse_level (l) {
      {foreach_child()
        val(w,0,0,0) = val(s,0,0,0);end_foreach_child()}
      _attribute[s.i].prolongation (point, s);
      {foreach_child() {
        double sp = val(s,0,0,0);
        val(s,0,0,0) = val(w,0,0,0);

        val(w,0,0,0) -= sp;
      }end_foreach_child()}
    }end_foreach_coarse_level();}END_FOREACH 
    boundary_level (((scalar[]){w,{-1}}), l + 1);
  }

  foreach_level_stencil(1,DEPAREN({0}),)
    { _stencil_val(s,0,0,0);_stencil_val_a(w,0,0,0); }end_foreach_level_stencil()

   BEGIN_FOREACH{
#line 104
foreach_level(0)
    val(w,0,0,0) = val(s,0,0,0);end_foreach_level();}END_FOREACH 
  boundary_level (((scalar[]){w,{-1}}), 0);
}

void inverse_wavelet (scalar s, scalar w)
{
  foreach_level_stencil(1,DEPAREN({0}),)
    { _stencil_val(w,0,0,0);_stencil_val_a(s,0,0,0); }end_foreach_level_stencil()
   BEGIN_FOREACH{
#line 111
foreach_level(0)
    val(s,0,0,0) = val(w,0,0,0);end_foreach_level();}END_FOREACH 
  boundary_level (((scalar[]){s,{-1}}), 0);
  for (int l = 0; l <= grid->maxdepth - 1; l++) {
    foreach_coarse_level_stencil (1,DEPAREN({l}),) {
default_stencil (
      
#line 116
point,((scalar[]){ s,{-1}}));
      {foreach_child()
        { _stencil_val(w,0,0,0);_stencil_val_r(s,0,0,0); }end_foreach_child()}
    }end_foreach_coarse_level_stencil()
     BEGIN_FOREACH{
#line 115
foreach_coarse_level (l) {
      _attribute[s.i].prolongation (point, s);
      {foreach_child()
        val(s,0,0,0) += val(w,0,0,0);end_foreach_child()}
    }end_foreach_coarse_level();}END_FOREACH 
    boundary_level (((scalar[]){s,{-1}}), l + 1);
  }
}

static inline double bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;







    return (27.*coarse(s,0,0,0) +
     9.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0) +
  coarse(s,0,0,child.z)) +
     3.*(coarse(s,child.x,child.y,0) + coarse(s,child.x,0,child.z) +
  coarse(s,0,child.y,child.z)) +
     coarse(s,child.x,child.y,child.z))/64.;

}


#line 124
static void _stencil_bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;







_stencil_coarse(s,0,0,0);
_stencil_coarse(s,o_stencil,0,0); _stencil_coarse(s,0,o_stencil,0);
  _stencil_coarse(s,0,0,o_stencil);
_stencil_coarse(s,o_stencil,o_stencil,0); _stencil_coarse(s,o_stencil,0,o_stencil);
  _stencil_coarse(s,0,o_stencil,o_stencil);
     _stencil_coarse(s,o_stencil,o_stencil,o_stencil);







    
#line 133
return 
        
         


;

}

static inline void refine_bilinear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = bilinear (point, s);end_foreach_child()}
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
#line 169 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid-common.h"
  if (!(false)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/multigrid-common.h", 169, "false");
  return 0.;

}

static inline double biquadratic_vertex (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;






  if (!(false)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/multigrid-common.h", 182, "false");
  return 0.;

}

static inline void refine_biquadratic (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  {foreach_child()
    val(s,0,0,0) = biquadratic (point, s);end_foreach_child()}
}

static inline void refine_linear (Point point, scalar s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;

#line 194
if(!is_constant(cm)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
      
#line 198
g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
      
#line 201
g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val(cm,0,0,0), sum = val(cm,0,0,0)*(1 << 3);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*val(cm,-child.x,0,0)/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*val(cm,0,-child.y,0)/cmc;
      
#line 207
val(s,0,0,0) += child.z*g.z*val(cm,0,0,-child.z)/cmc;
    sum -= val(cm,0,0,0);
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
}}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

#line 194
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  coord g;
  if (_attribute[s.i].gradient)
    {
      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
      
#line 198
g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
      
#line 198
g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
      
#line 201
g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
      
#line 201
g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*_const_cm, sum = _const_cm*(1 << 3);
  {foreach_child() {
    val(s,0,0,0) = sc;
    
      val(s,0,0,0) += child.x*g.x*_const_cm/cmc;
      
#line 207
val(s,0,0,0) += child.y*g.y*_const_cm/cmc;
      
#line 207
val(s,0,0,0) += child.z*g.z*_const_cm/cmc;
    sum -= _const_cm;
  }end_foreach_child()}
  if (!(fabs(sum) < 1e-10)) qassert ("/home/esteban/ProgramFile/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
}}

#line 211
}

static inline void refine_reset (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  {foreach_child()
    val(v,0,0,0) = 0.;end_foreach_child()}
}

static inline void refine_injection (Point point, scalar v)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double val = val(v,0,0,0);
  {foreach_child()
    val(v,0,0,0) = val;end_foreach_child()}
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static scalar multigrid_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_vertex_scalar (s, name);
  _attribute[s.i].restriction = restriction_vertex;
  return s;
}

static void multigrid_setup_vector (vector v)
{
   {
    _attribute[v.x.i].prolongation = refine_bilinear;
    _attribute[v.x.i].restriction = restriction_average;
  } 
#line 243
{
    _attribute[v.y.i].prolongation = refine_bilinear;
    _attribute[v.y.i].restriction = restriction_average;
  } 
#line 243
{
    _attribute[v.z.i].prolongation = refine_bilinear;
    _attribute[v.z.i].restriction = restriction_average;
  }
}

static vector multigrid_init_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  multigrid_setup_vector (v);
  return v;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  
    _attribute[v.y.i].restriction = no_restriction;
    
#line 260
_attribute[v.z.i].restriction = no_restriction;
    
#line 260
_attribute[v.x.i].restriction = no_restriction;
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

static tensor multigrid_init_tensor (tensor t, const char * name)
{
  t = cartesian_init_tensor (t, name);
  
    multigrid_setup_vector (t.x);
    
#line 269
multigrid_setup_vector (t.y);
    
#line 269
multigrid_setup_vector (t.z);
  return t;
}

void multigrid_debug (Point point)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 307 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      double zc = z - child.z*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++)
   for (int m = 0; m <= 1; m++) {
     {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){
       fprintf (fp, "%g %g %g %g ",
         xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
         yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
         zc + m*child.z*Delta*2. + _attribute[v.i].d.z*Delta,
         coarse(v,k*child.x,l*child.y,m*child.z));}}
     fputc ('\n', fp);
   }
      fprintf (ferr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);

    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 364 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4., zf = z - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++)
   for (int m = -2; m <= 3; m++) {
     {scalar*_i=(scalar*)( all);if(_i)for(scalar v=*_i;(&v)->i>=0;v=*++_i){ {
       fprintf (fp, "%g %g %g ",
         xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
         yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.,
         zf + m*Delta/2. + _attribute[v.i].d.z*Delta/4.);
       if (allocated_child(k,l,m))
  fprintf (fp, "%g ", fine(v,k,l,m));
       else
  fputs ("n/a ", fp);
     }}}
     fputc ('\n', fp);
   }
      fprintf (ferr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);

    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant (s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
     list2 = list_add (list2, _attribute[s.i].v.x);
     
#line 404
list2 = list_add (list2, _attribute[s.i].v.y);
     
#line 404
list2 = list_add (list2, _attribute[s.i].v.z);}
 else
   list2 = list_add (list2, s);
      }
    }}}

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
      foreach_coarse_level_stencil(1,DEPAREN({l}),) {
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     _stencil_restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {

default_stencil (
  
     
#line 418
point,((scalar[]){ s,{-1}}));
 }}}
      }end_foreach_coarse_level_stencil()
       BEGIN_FOREACH{
#line 412
foreach_coarse_level(l) {
 {scalar*_i=(scalar*)( listdef);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     restriction_average (point, s);}}
 {scalar*_i=(scalar*)( listc);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
  
     _attribute[s.i].restriction (point, s);
 }}}
      }end_foreach_coarse_level();}END_FOREACH 
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  init_scalar = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_vector = multigrid_init_vector;
  init_face_vector = multigrid_init_face_vector;
  init_tensor = multigrid_init_tensor;
  restriction = multigrid_restriction;
  debug = multigrid_debug;
}







void subtree_size (scalar size, bool leaves)
{




  foreach_stencil(1,{0},)
    {_stencil_val_a(size,0,0,0);  }end_foreach_stencil()




   BEGIN_FOREACH{
#line 453
foreach()
    val(size,0,0,0) = 1;end_foreach();}END_FOREACH 





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
    foreach_coarse_level_stencil(1,DEPAREN({l}),) {   
      
      {foreach_child()
 { _stencil_val(size,0,0,0); }end_foreach_child()}
      _stencil_val_a(size,0,0,0);  
    }end_foreach_coarse_level_stencil()
     BEGIN_FOREACH{
#line 462
foreach_coarse_level(l) {
      double sum = !leaves;
      {foreach_child()
 sum += val(size,0,0,0);end_foreach_child()}
      val(size,0,0,0) = sum;
    }end_foreach_coarse_level();}END_FOREACH 
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b,((scalar[]) {size,{-1}}), l); };
  }
}
#line 863 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"


void dimensions (int nx, int ny, int nz)
{

  int p[] = {nx, ny, nz};
  for (int i = 0; i < 3; i++)
    mpi_dims[i] = p[i];

}



#if 3 == 1

#define foreach_slice_x(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = start; point.i < end; point.i++)\

#line 882

#define end_foreach_slice_x() }

#elif 3 == 2

#define foreach_slice_x(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = start; point.i < end; point.i++)\
    for (point.j = 0; point.j < point.n + 2*2; point.j++)\

#line 892

#define end_foreach_slice_x() }

#define foreach_slice_y(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = 0; point.i < point.n + 2*2; point.i++)\
    for (point.j = start; point.j < end; point.j++)\

#line 900

#define end_foreach_slice_y() }

#elif 3 == 3

#define foreach_slice_x(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = start; point.i < end; point.i++)\
    for (point.j = 0; point.j < point.n + 2*2; point.j++)\
      for (point.k = 0; point.k < point.n + 2*2; point.k++)\

#line 911

#define end_foreach_slice_x() }

#define foreach_slice_y(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = 0; point.i < point.n + 2*2; point.i++)\
    for (point.j = start; point.j < end; point.j++)\
      for (point.k = 0; point.k < point.n + 2*2; point.k++)\

#line 920

#define end_foreach_slice_y() }

#define foreach_slice_z(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = 0; point.i < point.n + 2*2; point.i++)\
    for (point.j = 0; point.j < point.n + 2*2; point.j++)\
      for (point.k = start; point.k < end; point.k++)\

#line 929

#define end_foreach_slice_z() }

#endif

#line 1 "grid/multigrid-mpi.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid-mpi.h"


typedef struct {
  Boundary b;
  MPI_Comm cartcomm;
} MpiBoundary;


static void * snd_x (int i, int dst, int tag, int level, scalar * list,
       MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
   BEGIN_FOREACH{foreach_slice_x (i, i + 2, level)
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    }}}end_foreach_slice_x();}END_FOREACH 
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}

#line 9
static void * snd_y (int i, int dst, int tag, int level, scalar * list,
       MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
   BEGIN_FOREACH{foreach_slice_y (i, i + 2, level)
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    }}}end_foreach_slice_y();}END_FOREACH 
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}

#line 9
static void * snd_z (int i, int dst, int tag, int level, scalar * list,
       MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
   BEGIN_FOREACH{foreach_slice_z (i, i + 2, level)
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    }}}end_foreach_slice_z();}END_FOREACH 
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}


static void rcv_x (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s);
   BEGIN_FOREACH{foreach_slice_x (i, i + 2, level)
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    }}}end_foreach_slice_x();}END_FOREACH 
  pfree (buf,__func__,__FILE__,__LINE__);
}

#line 29
static void rcv_y (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s);
   BEGIN_FOREACH{foreach_slice_y (i, i + 2, level)
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    }}}end_foreach_slice_y();}END_FOREACH 
  pfree (buf,__func__,__FILE__,__LINE__);
}

#line 29
static void rcv_z (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    size += _attribute[s.i].block;}}
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s);
   BEGIN_FOREACH{foreach_slice_z (i, i + 2, level)
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    }}}end_foreach_slice_z();}END_FOREACH 
  pfree (buf,__func__,__FILE__,__LINE__);
}

     
static void mpi_boundary_level (const Boundary * b, scalar * list, int level)
{tracing("mpi_boundary_level","/home/esteban/ProgramFile/basilisk/src/grid/multigrid-mpi.h",49);
  scalar * list1 = NULL;
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s) && _attribute[s.i].block > 0)
      list1 = list_add (list1, s);}}
  if (!list1)
    {end_tracing("mpi_boundary_level","/home/esteban/ProgramFile/basilisk/src/grid/multigrid-mpi.h",56);return;}

  prof_start ("mpi_boundary_level");

  if (level < 0) level = depth();
  MpiBoundary * mpi = (MpiBoundary *) b;
  struct { int x, y, z; } dir = {0,1,2};
   {
    int left, right;
    MPI_Cart_shift (mpi->cartcomm, dir.x, 1, &left, &right);
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*2, nr = 0;
    if ((buf[0] = snd_x (npl - 2*2, right, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_x (2, left, 1, level, list1, &reqs[nr])))
      nr++;
    rcv_x (0, left, 0, level, list1);
    rcv_x (npl - 2, right, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    pfree (buf[0],__func__,__FILE__,__LINE__); pfree (buf[1],__func__,__FILE__,__LINE__);
  } 
#line 63
{
    int bottom, top;
    MPI_Cart_shift (mpi->cartcomm, dir.y, 1, &bottom, &top);
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*2, nr = 0;
    if ((buf[0] = snd_y (npl - 2*2, top, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_y (2, bottom, 1, level, list1, &reqs[nr])))
      nr++;
    rcv_y (0, bottom, 0, level, list1);
    rcv_y (npl - 2, top, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    pfree (buf[0],__func__,__FILE__,__LINE__); pfree (buf[1],__func__,__FILE__,__LINE__);
  } 
#line 63
{
    int back, front;
    MPI_Cart_shift (mpi->cartcomm, dir.z, 1, &back, &front);
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*2, nr = 0;
    if ((buf[0] = snd_z (npl - 2*2, front, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_z (2, back, 1, level, list1, &reqs[nr])))
      nr++;
    rcv_z (0, back, 0, level, list1);
    rcv_z (npl - 2, front, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    pfree (buf[0],__func__,__FILE__,__LINE__); pfree (buf[1],__func__,__FILE__,__LINE__);
  }

  pfree (list1,__func__,__FILE__,__LINE__);

  prof_stop();
end_tracing("mpi_boundary_level","/home/esteban/ProgramFile/basilisk/src/grid/multigrid-mpi.h",83);}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  MPI_Comm_free (&m->cartcomm);
  pfree (m,__func__,__FILE__,__LINE__);
}

Boundary * mpi_boundary_new()
{
  MpiBoundary * m = ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,__LINE__));
  MPI_Dims_create (npe(), 3, mpi_dims);
  MPI_Cart_create (MPI_COMM_WORLD, 3,
     mpi_dims, &Period.x, 0, &m->cartcomm);
  MPI_Cart_coords (m->cartcomm, pid(), 3, mpi_coords);


  struct { int x, y, z; } dir = {0,1,2};
   {
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.x, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (left);
    if (r != MPI_PROC_NULL)
      periodic_boundary (right);
  } 
#line 102
{
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.y, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (bottom);
    if (r != MPI_PROC_NULL)
      periodic_boundary (top);
  } 
#line 102
{
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.z, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (back);
    if (r != MPI_PROC_NULL)
      periodic_boundary (front);
  }


  N /= mpi_dims[0];
  int r = 0;
  while (N > 1)
    N /= 2, r++;
  grid->depth = grid->maxdepth = r;
  N = mpi_dims[0]*(1 << r);
  grid->n = 1 << 3*depth();
  grid->tn = npe()*grid->n;


  Boundary * b = (Boundary *) m;
  b->level = mpi_boundary_level;
  b->destroy = mpi_boundary_destroy;
  add_boundary (b);

  return b;
}

     
double z_indexing (scalar index, bool leaves)
{tracing("z_indexing","/home/esteban/ProgramFile/basilisk/src/grid/multigrid-mpi.h",131);
  long i;
  if (leaves)
    i = pid()*(1 << 3*depth());
  else
    i = pid()*((1 << 3*(depth() + 1)) - 1)/((1 << 3) - 1);
   BEGIN_FOREACH{foreach_cell() {
    if (!leaves || is_leaf(cell))
      val(index,0,0,0) = i++;
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}END_FOREACH 
  boundary_internal ((scalar *)((scalar[]){index,{-1}}), "/home/esteban/ProgramFile/basilisk/src/grid/multigrid-mpi.h", 144);
  { double _ret= pid() == 0 ? i*npe() - 1 : -1;end_tracing("z_indexing","/home/esteban/ProgramFile/basilisk/src/grid/multigrid-mpi.h",145);return _ret;}
end_tracing("z_indexing","/home/esteban/ProgramFile/basilisk/src/grid/multigrid-mpi.h",146);}
#line 935 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid.h"
#line 3 "/home/esteban/ProgramFile/basilisk/src/grid/multigrid3D.h"

void multigrid3D_methods() {
  multigrid_methods();
}
#line 10 "main3D.c"

#line 1 "navier-stokes/centered.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
#line 27 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
#line 1 "./run.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/run.h"
#line 9 "/home/esteban/ProgramFile/basilisk/src/run.h"
double dt = 1.;

#line 1 "./utils.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/utils.h"







double DT = 1e30, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf = {0};





void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;

  s.avg = mpi_time - t.tm;

  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
    
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(+:n)){
#line 69
foreach() n++;end_foreach();mpi_all_reduce_array(&n,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
    
#line 70
s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;





  s.mem = 0;


  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();



  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{



  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Multigrid"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(real)));

  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);

  fflush (fout);
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
  if(!is_constant(cm)){
  
#line 139
foreach_stencil(1,{0},
)
    {_stencil_val(f,0,0,0);_stencil_val(cm,0,0,0); {   
      _stencil_val(f,0,0,0);   
         
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0); 
       
    
#line 147
}       }end_foreach_stencil()
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 139
foreach(
)
    if (val(f,0,0,0) != 1e30 && (cube(Delta)*val(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (cube(Delta)*val(cm,0,0,0));
      avg += (cube(Delta)*val(cm,0,0,0))*v;
      rms += (cube(Delta)*val(cm,0,0,0))*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 147
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 139
foreach_stencil(1,{0},
)
    {_stencil_val(f,0,0,0);; {   
      _stencil_val(f,0,0,0);

;
;
; 
       
    
#line 147
}       }end_foreach_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(+:volume)
   reduction(+:rms) reduction(+:avg)reduction(max:max)){
#line 139
foreach(
)
    if (val(f,0,0,0) != 1e30 && (cube(Delta)*_const_cm) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (cube(Delta)*_const_cm);
      avg += (cube(Delta)*_const_cm)*v;
      rms += (cube(Delta)*_const_cm)*sq(v);
    }end_foreach();mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&rms,double,MPI_SUM,1);mpi_all_reduce_array(&avg,double,MPI_SUM,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 147
}
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  if(!is_constant(cm)){
  
#line 167
foreach_stencil(1,{0},
)
    {_stencil_val(cm,0,0,0); _stencil_val(f,0,0,0); {
_stencil_val(cm,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0); 
       _stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
_stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
         
         
    
#line 175
}      }end_foreach_stencil()
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 167
foreach(
)
    if ((cube(Delta)*val(cm,0,0,0)) > 0. && val(f,0,0,0) != 1e30) {
      volume += (cube(Delta)*val(cm,0,0,0));
      sum += (cube(Delta)*val(cm,0,0,0))*val(f,0,0,0);
      sum2 += (cube(Delta)*val(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 175
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 167
foreach_stencil(1,{0},
)
    {; _stencil_val(f,0,0,0); {
;
;_stencil_val(f,0,0,0);
;_stencil_val(f,0,0,0); 
       _stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
_stencil_val(f,0,0,0); { _stencil_val(f,0,0,0); }
         
         
    
#line 175
}      }end_foreach_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(min:min)
   reduction(max:max) reduction(+:volume) reduction(+:sum2)reduction(+:sum)){
#line 167
foreach(
)
    if ((cube(Delta)*_const_cm) > 0. && val(f,0,0,0) != 1e30) {
      volume += (cube(Delta)*_const_cm);
      sum += (cube(Delta)*_const_cm)*val(f,0,0,0);
      sum2 += (cube(Delta)*_const_cm)*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    }end_foreach();mpi_all_reduce_array(&min,double,MPI_MIN,1);mpi_all_reduce_array(&max,double,MPI_MAX,1);mpi_all_reduce_array(&volume,double,MPI_SUM,1);mpi_all_reduce_array(&sum2,double,MPI_SUM,1);mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 175
}
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
#line 191 "/home/esteban/ProgramFile/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
#line 217 "/home/esteban/ProgramFile/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
#line 241 "/home/esteban/ProgramFile/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  if (!(list_len(f) == vectors_len(g))) qassert ("/home/esteban/ProgramFile/basilisk/src/utils.h", 243, "list_len(f) == vectors_len(g)");
  foreach_stencil(1,{0},) {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {





_stencil_val(s,-1,0,0); _stencil_val(s,0,0,0); _stencil_val(s,1,0,0);





     
#line 254
_stencil_val_a(v.x,0,0,0);   
 } 
#line 248
{





_stencil_val(s,0,-1,0); _stencil_val(s,0,0,0); _stencil_val(s,0,1,0);





     
#line 254
_stencil_val_a(v.y,0,0,0);   
 } 
#line 248
{





_stencil_val(s,0,0,-1); _stencil_val(s,0,0,0); _stencil_val(s,0,0,1);





     
#line 254
_stencil_val_a(v.z,0,0,0);   
 }}
      else
 { {





_stencil_val(s,1,0,0); _stencil_val(s,-1,0,0);





     
#line 263
_stencil_val_a(v.x,0,0,0);   
 } 
#line 257
{





_stencil_val(s,0,1,0); _stencil_val(s,0,-1,0);





     
#line 263
_stencil_val_a(v.y,0,0,0);   
 } 
#line 257
{





_stencil_val(s,0,0,1); _stencil_val(s,0,0,-1);





     
#line 263
_stencil_val_a(v.z,0,0,0);   
 }}
    }}}
  }end_foreach_stencil()
   BEGIN_FOREACH{
#line 244
foreach() {
    scalar s; vector v;
    {vector*_i0=g;scalar*_i1= f;if(_i0)for(v=*_i0,s=*_i1;_i0->x.i>= 0;v=*++_i0,s=*++_i1){ {
      if (_attribute[s.i].gradient)
 { {





     val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
 } 
#line 248
{





     val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
 } 
#line 248
{





     val(v.z,0,0,0) = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1))/Delta;
 }}
      else
 { {





     val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
 } 
#line 257
{





     val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
 } 
#line 257
{





     val(v.z,0,0,0) = (val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta);
 }}
    }}}
  }end_foreach();}END_FOREACH 
}
#line 284 "/home/esteban/ProgramFile/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{
  if(!is_constant(fm.x) && !is_constant(cm)){
  
#line 286
foreach_stencil(1,{0},)
    {_stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,0,0,0);
        _stencil_val(fm.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,-1,0,0);
_stencil_val(fm.y,0,1,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,0,0);
        _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,-1,0); _stencil_val(fm.y,0,1,0);_stencil_val(u.x,0,1,0);_stencil_val(cm,0,0,0);
#line 287
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil() BEGIN_FOREACH{
#line 286
foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 0.)*Delta);end_foreach();}END_FOREACH }else if(is_constant(fm.x) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 286
foreach_stencil(1,{0},)
    {;;_stencil_val(u.y,0,0,0);
;_stencil_val(u.y,1,0,0);;_stencil_val(u.y,-1,0,0);
;;_stencil_val(u.x,0,0,0);
;_stencil_val(u.x,0,-1,0);;_stencil_val(u.x,0,1,0);_stencil_val(cm,0,0,0);
#line 287
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil()
   BEGIN_FOREACH{
#line 286
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(val(cm,0,0,0) + 0.)*Delta);end_foreach();}END_FOREACH }else if(!is_constant(fm.x) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 286
foreach_stencil(1,{0},)
    {_stencil_val(fm.x,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,0,0,0);
        _stencil_val(fm.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(fm.x,0,0,0);_stencil_val(u.y,-1,0,0);
_stencil_val(fm.y,0,1,0); _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,0,0);
        _stencil_val(fm.y,0,0,0);_stencil_val(u.x,0,-1,0); _stencil_val(fm.y,0,1,0);_stencil_val(u.x,0,1,0);;
#line 287
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil()
   BEGIN_FOREACH{
#line 286
foreach()
    val(omega,0,0,0) = ((val(fm.x,1,0,0) - val(fm.x,0,0,0))*val(u.y,0,0,0) +
        val(fm.x,1,0,0)*val(u.y,1,0,0) - val(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val(fm.y,0,1,0) - val(fm.y,0,0,0))*val(u.x,0,0,0) +
        val(fm.y,0,0,0)*val(u.x,0,-1,0) - val(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*(_const_cm + 0.)*Delta);end_foreach();}END_FOREACH }else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 286
foreach_stencil(1,{0},)
    {;;_stencil_val(u.y,0,0,0);
;_stencil_val(u.y,1,0,0);;_stencil_val(u.y,-1,0,0);
;;_stencil_val(u.x,0,0,0);
;_stencil_val(u.x,0,-1,0);;_stencil_val(u.x,0,1,0);;
#line 287
_stencil_val_a(omega,0,0,0);      
             

}end_foreach_stencil()
   BEGIN_FOREACH{
#line 286
foreach()
    val(omega,0,0,0) = ((_const_fm.x - _const_fm.x)*val(u.y,0,0,0) +
        _const_fm.x*val(u.y,1,0,0) - _const_fm.x*val(u.y,-1,0,0) -
        (_const_fm.y - _const_fm.y)*val(u.x,0,0,0) +
        _const_fm.y*val(u.x,0,-1,0) - _const_fm.y*val(u.x,0,1,0))/(2.*(_const_cm + 0.)*Delta);end_foreach();}END_FOREACH }
}





double change (scalar s, scalar sn)
{
  double max = 0.;
  if(!is_constant(cm)){
  
#line 300
foreach_stencil(1,{0},) {
_stencil_val(cm,0,0,0); {     
       _stencil_val(sn,0,0,0);_stencil_val(s,0,0,0);   
       
  
    } 
_stencil_val(s,0,0,0);
       
    
#line 306
_stencil_val_a(sn,0,0,0); 
  }end_foreach_stencil()
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(max:max)){
#line 300
foreach() {
    if ((cube(Delta)*val(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 307
}else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 300
foreach_stencil(1,{0},) {
; {     
       _stencil_val(sn,0,0,0);_stencil_val(s,0,0,0);   
       
  
    } 
_stencil_val(s,0,0,0);
       
    
#line 306
_stencil_val_a(sn,0,0,0); 
  }end_foreach_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(max:max)){
#line 300
foreach() {
    if ((cube(Delta)*_const_cm) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  }end_foreach();mpi_all_reduce_array(&max,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 307
}
  return max;
}





scalar lookup_field (const char * name)
{
  if (name)
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, name))
 return s;}}
  return (scalar){-1};
}

vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!strcmp (_attribute[s.i].name, component))
 return _attribute[s.i].v;}}
  }
  return (vector){{-1}};
}
#line 344 "/home/esteban/ProgramFile/basilisk/src/utils.h"
#define foreach_segment(_S,_p) {\
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};\
  double norm = sqrt(sq(t.x) + sq(t.y));\
  if (!(norm > 0.)) qassert ("/home/esteban/ProgramFile/basilisk/src/utils.h", 347, "norm > 0.");\
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;\
  double alpha = ((_S)[0].x*((_S)[1].y - (_S)[0].y) -\
    (_S)[0].y*((_S)[1].x - (_S)[0].x))/norm;\
  foreach()\
    if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {\
      coord _o = {x,y}, _p[2];\
      int _n = 0;\
 if (t.x)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].x = _o.x + _i*Delta/2.;\
     double a = (_p[_n].x - (_S)[0].x)/t.x;\
     _p[_n].y = (_S)[0].y + a*t.y;\
     if (fabs(_p[_n].y - _o.y) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].x = (_S)[0].x + a*t.x, _p[_n].y = (_S)[0].y + a*t.y;\
       if (fabs(_p[_n].x - _o.x) <= Delta/2. &&\
    fabs(_p[_n].y - _o.y) <= Delta/2.)\
  _n++;\
     }\
   }\
\
 if (t.y)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].y = _o.y + _i*Delta/2.;\
     double a = (_p[_n].y - (_S)[0].y)/t.y;\
     _p[_n].x = (_S)[0].x + a*t.x;\
     if (fabs(_p[_n].x - _o.x) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].y = (_S)[0].y + a*t.y, _p[_n].x = (_S)[0].x + a*t.x;\
       if (fabs(_p[_n].y - _o.y) <= Delta/2. &&\
    fabs(_p[_n].x - _o.x) <= Delta/2.)\
  _n++;\
     }\
   }\
\
      if (_n == 2) {\

#line 384

#line 414 "/home/esteban/ProgramFile/basilisk/src/utils.h"
#define end_foreach_segment() } } end_foreach(); }




void fields_stats()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    fprintf (ferr, " %s", _attribute[s.i].name);}}
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
    "name", "min", "avg", "stddev", "max");
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
      _attribute[s.i].name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }}}
}

#line 1 "./output.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/output.h"
#line 37 "/home/esteban/ProgramFile/basilisk/src/output.h"
     
void output_field (scalar * list,
     FILE * fp,
     int n,
     bool linear,
     coord box[2])
{tracing("output_field","/home/esteban/ProgramFile/basilisk/src/output.h",38);
  n++;
  int len = list_len (list);
  double Delta = 0.999999*(box[1].x - box[0].x)/(n - 1);
  int ny = (box[1].y - box[0].y)/Delta + 1;
  double ** field = (double **) matrix_new (n, ny, len*sizeof(double)), * v = field[0];
  for (int i = 0; i < n*ny*len; i++, v++)
    *v = 1e30;
  coord box1[2] = {{box[0].x - Delta/2., box[0].y - Delta/2.},
     {box[0].x + (n - 0.5)*Delta, box[0].y + (ny - 0.5)*Delta}};
  coord cn = {n, ny}, p;

  v = field[0];
  foreach_region_stencil (1,DEPAREN({p, box1, cn}),)



  {                     
    
    
    
    
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      { _stencil_interpolate_linear (point, s, NULL, NULL, NULL); _stencil_val(s,0,0,0);      }}}
  }end_foreach_region_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(min:v[:n*ny*len])){
#line 56
foreach_region (p, box1, cn)



  {
    double ** alias = field;
    int i = (p.x - box1[0].x)/(box1[1].x - box1[0].x)*cn.x;
    int j = (p.y - box1[0].y)/(box1[1].y - box1[0].y)*cn.y;
    int k = 0;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      alias[i][len*j + k++] = linear ? interpolate_linear (point, s, p.x, p.y, p.z) : val(s,0,0,0);}}
  }end_foreach_region();mpi_all_reduce_array(v,double,MPI_MIN,n*ny*len);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 

  
#line 69
if (pid() == 0) {
    fprintf (fp, "# 1:x 2:y");
    int i = 3;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      fprintf (fp, " %d:%s", i++, _attribute[s.i].name);}}
    fputc('\n', fp);
    for (int i = 0; i < n; i++) {
      double x = Delta*i + box[0].x;
      for (int j = 0; j < ny; j++) {
 double y = Delta*j + box[0].y;

 fprintf (fp, "%g %g", x, y);
 int k = 0;
 {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   fprintf (fp, " %g", field[i][len*j + k++]);}}
 fputc ('\n', fp);
      }
      fputc ('\n', fp);
    }
    fflush (fp);
  }

  matrix_free (field);
end_tracing("output_field","/home/esteban/ProgramFile/basilisk/src/output.h",92);}
#line 120 "/home/esteban/ProgramFile/basilisk/src/output.h"
     
void output_matrix (scalar f, FILE * fp, int n, bool linear)
{tracing("output_matrix","/home/esteban/ProgramFile/basilisk/src/output.h",121);
  float fn = n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, fp);
  for (int j = 0; j < n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, fp);
  }
  for (int i = 0; i < n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, fp);
    for (int j = 0; j < n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      v = interpolate (f, xp, yp
#line 899 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
, 0.
#line 135 "/home/esteban/ProgramFile/basilisk/src/output.h"
, linear);
      fwrite (&v, sizeof(float), 1, fp);
    }
  }
  fflush (fp);
end_tracing("output_matrix","/home/esteban/ProgramFile/basilisk/src/output.h",140);}
#line 149 "/home/esteban/ProgramFile/basilisk/src/output.h"
typedef void (* Colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

void blue_white_red (double cmap[127][3])
{
  for (int i = 0; i < (127 + 1)/2; i++) {
    cmap[i][0] = i/((127 - 1)/2.);
    cmap[i][1] = i/((127 - 1)/2.);
    cmap[i][2] = 1.;
  }
  for (int i = 0; i < (127 - 1)/2; i++) {
    cmap[i + (127 + 1)/2][0] = 1.;
    cmap[i + (127 + 1)/2][1] = cmap[(127 - 3)/2 - i][1];
    cmap[i + (127 + 1)/2][2] = cmap[(127 - 3)/2 - i][1];
  }
}





typedef struct {
  unsigned char r, g, b;
} Color;

Color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  Color c;
  if (val == 1e30) {
    c.r = c.g = c.b = 0;
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = 127 - 2, coef = 1.;
  else {
    i = val*(127 - 1);
    coef = val*(127 - 1) - i;
  }
  if (!(i >= 0 && i < 127 - 1)) qassert ("/home/esteban/ProgramFile/basilisk/src/output.h", 281, "i >= 0 && i < NCMAP - 1");
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 300 "/home/esteban/ProgramFile/basilisk/src/output.h"
static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup()
{
  for (int i = 0; i < open_image_data.n; i++) {
    qpclose (open_image_data.fp[i]);
    pfree (open_image_data.names[i],__func__,__FILE__,__LINE__);
  }
  pfree (open_image_data.fp,__func__,__FILE__,__LINE__);
  pfree (open_image_data.names,__func__,__FILE__,__LINE__);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);

    MPI_Abort (MPI_COMM_WORLD, 1);

    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  if (!(pid() == 0)) qassert ("/home/esteban/ProgramFile/basilisk/src/output.h", 382, "pid() == 0");
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
 has_ffmpeg = true;
      else {
 fprintf (ferr,
   "src/output.h:%d: warning: cannot find '%s' or 'ffmpeg'/'avconv'\n"
   "src/output.h:%d: warning: falling back to raw PPM outputs\n",
   402, command, 402);
 has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }
    open_image_data.n++;
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,__LINE__);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,__LINE__);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    return open_image_data.fp[open_image_data.n - 1] = qpopen (command, "w");
  }
  else {
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
 has_convert = true;
      else {
 fprintf (ferr,
   "src/output.h:%d: warning: cannot find 'convert'\n"
   "src/output.h:%d: warning: falling back to raw PPM outputs\n",
   436, 436);
 has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");

    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return qpopen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  if (!(pid() == 0)) qassert ("/home/esteban/ProgramFile/basilisk/src/output.h", 458, "pid() == 0");
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
#line 539 "/home/esteban/ProgramFile/basilisk/src/output.h"
     
void output_ppm (scalar f,
   FILE * fp,
   int n,
   char * file,
   double min, double max, double spread,
   double z,
   bool linear,
   coord box[2],
   scalar mask,
   Colormap map,
   char * opt,
   int fps)
{tracing("output_ppm","/home/esteban/ProgramFile/basilisk/src/output.h",540);

  if (!min && !max) {
    stats s = statsf (f);
    if (spread < 0.)
      min = s.min, max = s.max;
    else {
      double avg = s.sum/s.volume;
      min = avg - spread*s.stddev; max = avg + spread*s.stddev;
    }
  }
  box[0].z = z, box[1].z = z;

  coord cn = {n}, p;
  double delta = (box[1].x - box[0].x)/n;
  cn.y = (int)((box[1].y - box[0].y)/delta);
  if (((int)cn.y) % 2) cn.y++;

  Color ** ppm = (Color **) matrix_new (cn.y, cn.x, sizeof(Color));
  unsigned char * ppm0 = &ppm[0][0].r;
  int len = 3*cn.x*cn.y;
  memset (ppm0, 0, len*sizeof (unsigned char));
  double cmap[127][3];
  (* map) (cmap);


  foreach_region_stencil (1,DEPAREN({p, box, cn}),)



  { 
    
    if (mask.i >= 0) {
      if (linear) {  
  _stencil_interpolate_linear (point, mask, NULL, NULL, NULL);
{ 
    
   
{ _stencil_interpolate_linear (point, f, NULL, NULL, NULL); }}
    
 
      
#line 591
}
      else {
_stencil_val(mask,0,0,0);{
     
   
{ _stencil_val(f,0,0,0); }}
    
 
      
#line 597
}
    }
    else if (linear)
      { _stencil_interpolate_linear (point, f, NULL, NULL, NULL); }
    else
      { _stencil_val(f,0,0,0); }                  
    
    
         
         
  }end_foreach_region_stencil()


  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(max:ppm0[:len])){
#line 578
foreach_region (p, box, cn)



  {
    double v;
    if (mask.i >= 0) {
      if (linear) {
 double m = interpolate_linear (point, mask, p.x, p.y, p.z);
 if (m < 0.)
   v = 1e30;
 else
   v = interpolate_linear (point, f, p.x, p.y, p.z);
      }
      else {
 if (val(mask,0,0,0) < 0.)
   v = 1e30;
 else
   v = val(f,0,0,0);
      }
    }
    else if (linear)
      v = interpolate_linear (point, f, p.x, p.y, p.z);
    else
      v = val(f,0,0,0);
    int i = (p.x - box[0].x)/(box[1].x - box[0].x)*cn.x;
    int j = (p.y - box[0].y)/(box[1].y - box[0].y)*cn.y;
    Color ** alias = ppm;
    alias[(int)cn.y - 1 - j][i] = colormap_color (cmap, v, min, max);
  }end_foreach_region();mpi_all_reduce_array(ppm0,unsigned char,MPI_MAX,len);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 

  
#line 609
if (pid() == 0) {
    if (file)
      fp = open_image (file, opt);

    fprintf (fp, "P6\n%g %g 255\n", cn.x, cn.y);
    fwrite (ppm0, sizeof(unsigned char), 3*cn.x*cn.y, fp);

    if (file)
      close_image (file, fp);
    else
      fflush (fp);
  }

  matrix_free (ppm);
end_tracing("output_ppm","/home/esteban/ProgramFile/basilisk/src/output.h",623);}
#line 655 "/home/esteban/ProgramFile/basilisk/src/output.h"
     
void output_grd (scalar f,
   FILE * fp,
   double Delta,
   bool linear,
   double box[2][2],
   scalar mask)
{tracing("output_grd","/home/esteban/ProgramFile/basilisk/src/output.h",656);
  int nx = (box[1][0] - box[0][0])/Delta;
  int ny = (box[1][1] - box[0][1])/Delta;


  fprintf (fp, "ncols          %d\n", nx);
  fprintf (fp, "nrows          %d\n", ny);
  fprintf (fp, "xllcorner      %g\n", box[0][0]);
  fprintf (fp, "yllcorner      %g\n", box[0][1]);
  fprintf (fp, "cellsize       %g\n", Delta);
  fprintf (fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + box[0][0] + Delta/2., v;
      if (mask.i >= 0) {
 double m = interpolate (mask, xp, yp
#line 899 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
, 0.
#line 680 "/home/esteban/ProgramFile/basilisk/src/output.h"
, linear);
 if (m < 0.)
   v = 1e30;
 else
   v = interpolate (f, xp, yp
#line 899 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
, 0.
#line 684 "/home/esteban/ProgramFile/basilisk/src/output.h"
, linear);
      }
      else
 v = interpolate (f, xp, yp
#line 899 "/home/esteban/ProgramFile/basilisk/src/grid/cartesian-common.h"
, 0.
#line 687 "/home/esteban/ProgramFile/basilisk/src/output.h"
, linear);
      if (v == 1e30)
 fprintf (fp, "-9999 ");
      else
 fprintf (fp, "%f ", v);
    }
    fprintf (fp, "\n");
  }

  fflush (fp);
end_tracing("output_grd","/home/esteban/ProgramFile/basilisk/src/output.h",697);}
#line 724 "/home/esteban/ProgramFile/basilisk/src/output.h"
static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}

     
void output_gfs (FILE * fp,
   scalar * list,
   char * file,
   bool translate)
{tracing("output_gfs","/home/esteban/ProgramFile/basilisk/src/output.h",745);
  char * fname = file;



  not_mpi_compatible();

  FILE * sfp = fp;
  if (file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,__LINE__));
    snprintf (fname, 80, ".output-%ld", pid);
    fp = NULL;
  }


  bool opened = false;
  if (fp == NULL) {
    if (fname == NULL)
      fp = fout;
    else if (!(fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * slist = list ? list : list_copy (all);

  restriction (slist);
  fprintf (fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);

  fprintf (fp, "z = %g ", 0.5 + Z0/L0);


  if (slist != NULL && slist[0].i != -1) {
    scalar s = slist[0];
    char * name = replace (_attribute[s.i].name, '.', '_', translate);
    fprintf (fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(slist); i++) {
      scalar s = slist[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', translate);
 fprintf (fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (fp, " ");
  }
  fprintf (fp, "} {\n");
  fprintf (fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (fp, "  VariableTracerVOF f\n");
  fprintf (fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");


  long header;
  if ((header = ftell (fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (_attribute[s.i].name)
      cell_size += sizeof(double);}}
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;




   BEGIN_FOREACH{foreach_cell() {

    if (is_local(cell))

    {

      if (fseek (fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }

      unsigned flags =
 level == 0 ? 0 :
#line 848 "/home/esteban/ProgramFile/basilisk/src/output.h"
      child.x == -1 && child.y == -1 && child.z == -1 ? 0 :
 child.x == -1 && child.y == -1 && child.z == 1 ? 1 :
 child.x == -1 && child.y == 1 && child.z == -1 ? 2 :
 child.x == -1 && child.y == 1 && child.z == 1 ? 3 :
 child.x == 1 && child.y == -1 && child.z == -1 ? 4 :
 child.x == 1 && child.y == -1 && child.z == 1 ? 5 :
 child.x == 1 && child.y == 1 && child.z == -1 ? 6 :
 7;

      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, fp);
      {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != 1e30 ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != 1e30 ? - val(s,0,0,0) : (double) DBL_MAX;
     }


     else
       a = is_local(cell) && val(s,0,0,0) != 1e30 ? val(s,0,0,0) : (double) DBL_MAX;

   }
   else
     a = is_local(cell) && val(s,0,0,0) != 1e30 ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, fp);
 }}}
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}END_FOREACH 


  delete (((scalar[]){index,{-1}}));
  if (!pid() && fseek (fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())

    fputs ("}\n", fp);
  fflush (fp);

  if (!list)
    pfree (slist,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (fp);


  if (file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (sfp == NULL)
 sfp = fout;
      fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, fp)) > 0)
 fwrite (buffer, 1, l, sfp);
      fflush (sfp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,__LINE__);
  }

end_tracing("output_gfs","/home/esteban/ProgramFile/basilisk/src/output.h",925);}
#line 949 "/home/esteban/ProgramFile/basilisk/src/output.h"
struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar[]){cm,{-1}}), NULL);
  {scalar*_i=(scalar*)( lista);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);}}
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }}}
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}
#line 1046 "/home/esteban/ProgramFile/basilisk/src/output.h"
     
void dump (const char * file,
    scalar * list,
    FILE * fp,
    bool unbuffered)
{tracing("dump","/home/esteban/ProgramFile/basilisk/src/output.h",1047);
  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);
  }

  scalar * dlist = dump_list (list);
  scalar  size=new_scalar("size");
  scalar * slist = list_concat (((scalar[]){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(slist), iter, depth(), npe(),
          dump_version };


  for (int i = 0; i < 3; i++)
    (&header.n.x)[i] = mpi_dims[i];
  MPI_Barrier (MPI_COMM_WORLD);


  if (pid() == 0)
    dump_header (fh, &header, slist);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);}}
  long pos = pid() ? 0 : sizeofheader;

  subtree_size (size, false);

   BEGIN_FOREACH{foreach_cell() {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      if (pos != offset) {
 fseek (fh, offset, SEEK_SET);
 pos = offset;
      }
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
 double val = val(s,0,0,0);
 fwrite (&val, 1, sizeof(double), fh);
      }}}
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}END_FOREACH 

  delete (((scalar[]){index,{-1}}));

  pfree (slist,__func__,__FILE__,__LINE__);
  fclose (fh);
  if (!unbuffered && pid() == 0)
    rename (name, file);delete((scalar*)((scalar[]){size,{-1}}));
end_tracing("dump","/home/esteban/ProgramFile/basilisk/src/output.h",1120);}


     
bool restore (const char * file,
       scalar * list,
       FILE * fp)
{tracing("restore","/home/esteban/ProgramFile/basilisk/src/output.h",1124);
  if (!fp && (fp = fopen (file, "r")) == NULL)
    {end_tracing("restore","/home/esteban/ProgramFile/basilisk/src/output.h",1129);return false;}
  if (!(fp)) qassert ("/home/esteban/ProgramFile/basilisk/src/output.h", 1130, "fp");

  struct DumpHeader header = {0};
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }
#line 1147 "/home/esteban/ProgramFile/basilisk/src/output.h"
  if (header.npe != npe()) {
    fprintf (ferr,
      "restore(): error: the number of processes don't match:"
      " %d != %d\n",
      header.npe, npe());
    exit (1);
  }
  dimensions (header.n.x, header.n.y, header.n.z);
  double n = header.n.x;
  int depth = header.depth;
  while (n > 1)
    depth++, n /= 2;
  init_grid (1 << depth);





  bool restore_all = (list == all);
  scalar * slist = dump_list (list ? list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (slist)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (slist));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }}}
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (slist,__func__,__FILE__,__LINE__);
    slist = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin (o[0], o[1], o[2]);
    size (o[3]);
  }


  long cell_size = sizeof(unsigned) + header.len*sizeof(double);
  long offset = pid()*((1 << 3*(header.depth + 1)) - 1)/
    ((1 << 3) - 1)*cell_size;
  if (fseek (fp, offset, SEEK_CUR) < 0) {
    perror ("restore(): error while seeking");
    exit (1);
  }


  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector[]){fm,{{-1},{-1},{-1}}});



   BEGIN_FOREACH{foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }

    fseek (fp, sizeof(double), SEEK_CUR);
    {scalar*_i=(scalar*)( slist);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      double val;
      if (fread (&val, sizeof(double), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX)
 val(s,0,0,0) = val;
    }}}
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, listm, 0, NULL);
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}END_FOREACH 
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    _attribute[s.i].dirty = true;}}


  scalar * other = NULL;
  {scalar*_i=(scalar*)( all);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!list_lookup (slist, s) && !list_lookup (listm, s))
      other = list_append (other, s);}}
  reset (other, 0.);
  pfree (other,__func__,__FILE__,__LINE__);

  pfree (slist,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  {end_tracing("restore","/home/esteban/ProgramFile/basilisk/src/output.h",1290);return true;}
end_tracing("restore","/home/esteban/ProgramFile/basilisk/src/output.h",1291);}
#line 435 "/home/esteban/ProgramFile/basilisk/src/utils.h"
#line 12 "/home/esteban/ProgramFile/basilisk/src/run.h"

     
void run (void)
{tracing("run","/home/esteban/ProgramFile/basilisk/src/run.h",14);
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {





    update_perf();
    iter = inext, t = tnext;
  }




  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
end_tracing("run","/home/esteban/ProgramFile/basilisk/src/run.h",37);}




static int defaults_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}





#line 42
      static int defaults(const int i,const double t,Event *_ev){tracing("defaults","/home/esteban/ProgramFile/basilisk/src/run.h",42); {
  display ("box();"
#line 494 "/home/esteban/ProgramFile/basilisk/src/common.h"
, false
#line 43 "/home/esteban/ProgramFile/basilisk/src/run.h"
);
}{end_tracing("defaults","/home/esteban/ProgramFile/basilisk/src/run.h",44);return 0;}end_tracing("defaults","/home/esteban/ProgramFile/basilisk/src/run.h",44);}





static int cleanup_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = TEND_EVENT)!=0;*ip=i;*tp=t;return ret;}






#line 50
      static int cleanup(const int i,const double t,Event *_ev){tracing("cleanup","/home/esteban/ProgramFile/basilisk/src/run.h",50); {
  display ("", true);
}{end_tracing("cleanup","/home/esteban/ProgramFile/basilisk/src/run.h",52);return 0;}end_tracing("cleanup","/home/esteban/ProgramFile/basilisk/src/run.h",52);}
#line 28 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
#line 1 "./timestep.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  if (t == 0.) previous = 0.;
  dtmax /= CFL;
  if(!is_constant(fm.x)){
  
#line 7
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {_stencil_val(u.x,0,0,0); {   
      _stencil_val(u.x,0,0,0);
_stencil_val(fm.x,0,0,0); 
_stencil_val(fm.x,0,0,0);    
          
       
         
    
#line 13
}   }}end__stencil_is_face_x()
#line 7
_stencil_is_face_y(){
    {_stencil_val(u.y,0,0,0); {   
      _stencil_val(u.y,0,0,0);
_stencil_val(fm.y,0,0,0); 
_stencil_val(fm.y,0,0,0);    
          
       
         
    
#line 13
}   }}end__stencil_is_face_y()
#line 7
_stencil_is_face_z(){
    {_stencil_val(u.z,0,0,0); {   
      _stencil_val(u.z,0,0,0);
_stencil_val(fm.z,0,0,0); 
_stencil_val(fm.z,0,0,0);    
          
       
         
    
#line 13
}   }}end__stencil_is_face_z()}end_foreach_face_stencil()
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(min:dtmax)){
#line 7
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));
      if (!(val(fm.x,0,0,0))) qassert ("/home/esteban/ProgramFile/basilisk/src/timestep.h", 10, "fm.x[]");
      dt *= val(fm.x,0,0,0);
      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 7
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));
      if (!(val(fm.y,0,0,0))) qassert ("/home/esteban/ProgramFile/basilisk/src/timestep.h", 10, "fm.x[]");
      dt *= val(fm.y,0,0,0);
      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()
#line 7
is_face_z(){
    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.z,0,0,0));
      if (!(val(fm.z,0,0,0))) qassert ("/home/esteban/ProgramFile/basilisk/src/timestep.h", 10, "fm.x[]");
      dt *= val(fm.z,0,0,0);
      if (dt < dtmax) dtmax = dt;
    }}end_is_face_z()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 13
}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 7
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {_stencil_val(u.x,0,0,0); {   
      _stencil_val(u.x,0,0,0);
;
;    
          
       
         
    
#line 13
}   }}end__stencil_is_face_x()
#line 7
_stencil_is_face_y(){
    {_stencil_val(u.y,0,0,0); {   
      _stencil_val(u.y,0,0,0);
;
;    
          
       
         
    
#line 13
}   }}end__stencil_is_face_y()
#line 7
_stencil_is_face_z(){
    {_stencil_val(u.z,0,0,0); {   
      _stencil_val(u.z,0,0,0);
;
;    
          
       
         
    
#line 13
}   }}end__stencil_is_face_z()}end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(min:dtmax)){
#line 7
foreach_face_generic(){is_face_x(){
    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));
      if (!(_const_fm.x)) qassert ("/home/esteban/ProgramFile/basilisk/src/timestep.h", 10, "fm.x[]");
      dt *= _const_fm.x;
      if (dt < dtmax) dtmax = dt;
    }}end_is_face_x()
#line 7
is_face_y(){
    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));
      if (!(_const_fm.y)) qassert ("/home/esteban/ProgramFile/basilisk/src/timestep.h", 10, "fm.x[]");
      dt *= _const_fm.y;
      if (dt < dtmax) dtmax = dt;
    }}end_is_face_y()
#line 7
is_face_z(){
    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.z,0,0,0));
      if (!(_const_fm.z)) qassert ("/home/esteban/ProgramFile/basilisk/src/timestep.h", 10, "fm.x[]");
      dt *= _const_fm.z;
      if (dt < dtmax) dtmax = dt;
    }}end_is_face_z()}end_foreach_face_generic();mpi_all_reduce_array(&dtmax,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 13
}
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
#line 29 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
#line 1 "./bcg.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/bcg.h"
#line 11 "/home/esteban/ProgramFile/basilisk/src/bcg.h"
     
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
              scalar src)
{tracing("tracer_fluxes","/home/esteban/ProgramFile/basilisk/src/bcg.h",12);





  vector  g=new_vector("g");
  gradients (((scalar[]){f,{-1}}),((vector[]) {g,{{-1},{-1},{-1}}}));




  if(!is_constant(fm.x) && !is_constant(src)){




  
#line 29
foreach_face_stencil(1,{0},){_stencil_is_face_x(){ {        







    _stencil_val(fm.x,0,0,0);_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0); _stencil_val(src,-1,0,0);_stencil_val(src,0,0,0);_stencil_val(f, o_stencil,0,0);





_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {     
       _stencil_val(fm.y,o_stencil,1,0);_stencil_val(fm.y,o_stencil,0,0); _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }


_stencil_val(fm.z,o_stencil,0,0); _stencil_val(fm.z,o_stencil,0,1); {     
       _stencil_val(fm.z,o_stencil,0,1);_stencil_val(fm.z,o_stencil,0,0); _stencil_val(uf.z,o_stencil,0,1);_stencil_val(uf.z,o_stencil,0,0);         
       _stencil_val(f,o_stencil,0,-1);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,0,1);
        
    }


_stencil_val(uf.x,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 29
_stencil_is_face_y(){ {        







    _stencil_val(fm.y,0,0,0);_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0); _stencil_val(src,0,-1,0);_stencil_val(src,0,0,0);_stencil_val(f,0, o_stencil,0);





_stencil_val(fm.z,0,o_stencil,0); _stencil_val(fm.z,0,o_stencil,1); {     
       _stencil_val(fm.z,0,o_stencil,1);_stencil_val(fm.z,0,o_stencil,0); _stencil_val(uf.z,0,o_stencil,1);_stencil_val(uf.z,0,o_stencil,0);         
       _stencil_val(f,0,o_stencil,-1);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,0,o_stencil,1);
        
    }


_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {     
       _stencil_val(fm.x,1,o_stencil,0);_stencil_val(fm.x,0,o_stencil,0); _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }


_stencil_val(uf.y,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()
#line 29
_stencil_is_face_z(){ {        







    _stencil_val(fm.z,0,0,0);_stencil_val(uf.z,0,0,0);              
    
    _stencil_val(g.z,0,0,o_stencil); _stencil_val(src,0,0,-1);_stencil_val(src,0,0,0);_stencil_val(f,0,0, o_stencil);





_stencil_val(fm.x,0,0,o_stencil); _stencil_val(fm.x,1,0,o_stencil); {     
       _stencil_val(fm.x,1,0,o_stencil);_stencil_val(fm.x,0,0,o_stencil); _stencil_val(uf.x,1,0,o_stencil);_stencil_val(uf.x,0,0,o_stencil);         
       _stencil_val(f,-1,0,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,1,0,o_stencil);
        
    }


_stencil_val(fm.y,0,0,o_stencil); _stencil_val(fm.y,0,1,o_stencil); {     
       _stencil_val(fm.y,0,1,o_stencil);_stencil_val(fm.y,0,0,o_stencil); _stencil_val(uf.y,0,1,o_stencil);_stencil_val(uf.y,0,0,o_stencil);         
       _stencil_val(f,0,-1,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,0,1,o_stencil);
        
    }


_stencil_val(uf.z,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.z,0,0,0);  
  }}end__stencil_is_face_z()}end_foreach_face_stencil() BEGIN_FOREACH{
#line 29
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.z,i,0,0) && val(fm.z,i,0,1)) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(val(fm.z,i,0,0) + val(fm.z,i,0,1));
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 29
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.z,0,i,0) && val(fm.z,0,i,1)) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(val(fm.z,0,i,0) + val(fm.z,0,i,1));
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()
#line 29
is_face_z(){ {







    double un = dt*val(uf.z,0,0,0)/(val(fm.z,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val(src,0,0,0) + val(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (val(fm.x,0,0,i) && val(fm.x,1,0,i)) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(val(fm.x,0,0,i) + val(fm.x,1,0,i));
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.y,0,0,i) && val(fm.y,0,1,i)) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(val(fm.y,0,0,i) + val(fm.y,0,1,i));
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  }}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(is_constant(fm.x) && !is_constant(src)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);




  
#line 29
foreach_face_stencil(1,{0},){_stencil_is_face_x(){ {







;_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0); _stencil_val(src,-1,0,0);_stencil_val(src,0,0,0);_stencil_val(f, o_stencil,0,0);





;; {
;; _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }


;; {
;; _stencil_val(uf.z,o_stencil,0,1);_stencil_val(uf.z,o_stencil,0,0);         
       _stencil_val(f,o_stencil,0,-1);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,0,1);
        
    }


_stencil_val(uf.x,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 29
_stencil_is_face_y(){ {







;_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0); _stencil_val(src,0,-1,0);_stencil_val(src,0,0,0);_stencil_val(f,0, o_stencil,0);





;; {
;; _stencil_val(uf.z,0,o_stencil,1);_stencil_val(uf.z,0,o_stencil,0);         
       _stencil_val(f,0,o_stencil,-1);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,0,o_stencil,1);
        
    }


;; {
;; _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }


_stencil_val(uf.y,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()
#line 29
_stencil_is_face_z(){ {







;_stencil_val(uf.z,0,0,0);              
    
    _stencil_val(g.z,0,0,o_stencil); _stencil_val(src,0,0,-1);_stencil_val(src,0,0,0);_stencil_val(f,0,0, o_stencil);





;; {
;; _stencil_val(uf.x,1,0,o_stencil);_stencil_val(uf.x,0,0,o_stencil);         
       _stencil_val(f,-1,0,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,1,0,o_stencil);
        
    }


;; {
;; _stencil_val(uf.y,0,1,o_stencil);_stencil_val(uf.y,0,0,o_stencil);         
       _stencil_val(f,0,-1,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,0,1,o_stencil);
        
    }


_stencil_val(uf.z,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.z,0,0,0);  
  }}end__stencil_is_face_z()}end_foreach_face_stencil()




   BEGIN_FOREACH{
#line 29
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val(src,0,0,0) + val(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.z && _const_fm.z) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(_const_fm.z + _const_fm.z);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 29
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val(src,0,0,0) + val(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.z && _const_fm.z) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(_const_fm.z + _const_fm.z);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.x && _const_fm.x) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()
#line 29
is_face_z(){ {







    double un = dt*val(uf.z,0,0,0)/(_const_fm.z*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val(src,0,0,0) + val(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.y && _const_fm.y) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(_const_fm.y + _const_fm.y);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  }}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(!is_constant(fm.x) && is_constant(src)){double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  
#line 29
foreach_face_stencil(1,{0},){_stencil_is_face_x(){ {        







    _stencil_val(fm.x,0,0,0);_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0);;;_stencil_val(f, o_stencil,0,0);





_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {     
       _stencil_val(fm.y,o_stencil,1,0);_stencil_val(fm.y,o_stencil,0,0); _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }


_stencil_val(fm.z,o_stencil,0,0); _stencil_val(fm.z,o_stencil,0,1); {     
       _stencil_val(fm.z,o_stencil,0,1);_stencil_val(fm.z,o_stencil,0,0); _stencil_val(uf.z,o_stencil,0,1);_stencil_val(uf.z,o_stencil,0,0);         
       _stencil_val(f,o_stencil,0,-1);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,0,1);
        
    }


_stencil_val(uf.x,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 29
_stencil_is_face_y(){ {        







    _stencil_val(fm.y,0,0,0);_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0);;;_stencil_val(f,0, o_stencil,0);





_stencil_val(fm.z,0,o_stencil,0); _stencil_val(fm.z,0,o_stencil,1); {     
       _stencil_val(fm.z,0,o_stencil,1);_stencil_val(fm.z,0,o_stencil,0); _stencil_val(uf.z,0,o_stencil,1);_stencil_val(uf.z,0,o_stencil,0);         
       _stencil_val(f,0,o_stencil,-1);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,0,o_stencil,1);
        
    }


_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {     
       _stencil_val(fm.x,1,o_stencil,0);_stencil_val(fm.x,0,o_stencil,0); _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }


_stencil_val(uf.y,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()
#line 29
_stencil_is_face_z(){ {        







    _stencil_val(fm.z,0,0,0);_stencil_val(uf.z,0,0,0);              
    
    _stencil_val(g.z,0,0,o_stencil);;;_stencil_val(f,0,0, o_stencil);





_stencil_val(fm.x,0,0,o_stencil); _stencil_val(fm.x,1,0,o_stencil); {     
       _stencil_val(fm.x,1,0,o_stencil);_stencil_val(fm.x,0,0,o_stencil); _stencil_val(uf.x,1,0,o_stencil);_stencil_val(uf.x,0,0,o_stencil);         
       _stencil_val(f,-1,0,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,1,0,o_stencil);
        
    }


_stencil_val(fm.y,0,0,o_stencil); _stencil_val(fm.y,0,1,o_stencil); {     
       _stencil_val(fm.y,0,1,o_stencil);_stencil_val(fm.y,0,0,o_stencil); _stencil_val(uf.y,0,1,o_stencil);_stencil_val(uf.y,0,0,o_stencil);         
       _stencil_val(f,0,-1,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,0,1,o_stencil);
        
    }


_stencil_val(uf.z,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.z,0,0,0);  
  }}end__stencil_is_face_z()}end_foreach_face_stencil()




   BEGIN_FOREACH{
#line 29
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(val(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val(fm.y,i,0,0) + val(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.z,i,0,0) && val(fm.z,i,0,1)) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(val(fm.z,i,0,0) + val(fm.z,i,0,1));
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 29
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(val(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val(fm.z,0,i,0) && val(fm.z,0,i,1)) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(val(fm.z,0,i,0) + val(fm.z,0,i,1));
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val(fm.x,0,i,0) + val(fm.x,1,i,0));
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()
#line 29
is_face_z(){ {







    double un = dt*val(uf.z,0,0,0)/(val(fm.z,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (val(fm.x,0,0,i) && val(fm.x,1,0,i)) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(val(fm.x,0,0,i) + val(fm.x,1,0,i));
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val(fm.y,0,0,i) && val(fm.y,0,1,i)) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(val(fm.y,0,0,i) + val(fm.y,0,1,i));
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  }}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_src=_constant[src.i-_NVARMAX];NOT_UNUSED(_const_src);




  
#line 29
foreach_face_stencil(1,{0},){_stencil_is_face_x(){ {







;_stencil_val(uf.x,0,0,0);              
    
    _stencil_val(g.x,o_stencil,0,0);;;_stencil_val(f, o_stencil,0,0);





;; {
;; _stencil_val(uf.y,o_stencil,1,0);_stencil_val(uf.y,o_stencil,0,0);         
       _stencil_val(f,o_stencil,-1,0);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,1,0);
        
    }


;; {
;; _stencil_val(uf.z,o_stencil,0,1);_stencil_val(uf.z,o_stencil,0,0);         
       _stencil_val(f,o_stencil,0,-1);_stencil_val(f, o_stencil,0,0);_stencil_val(f, o_stencil,0,0); _stencil_val(f,o_stencil,0,1);
        
    }


_stencil_val(uf.x,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.x,0,0,0);  
  }}end__stencil_is_face_x()
#line 29
_stencil_is_face_y(){ {







;_stencil_val(uf.y,0,0,0);              
    
    _stencil_val(g.y,0,o_stencil,0);;;_stencil_val(f,0, o_stencil,0);





;; {
;; _stencil_val(uf.z,0,o_stencil,1);_stencil_val(uf.z,0,o_stencil,0);         
       _stencil_val(f,0,o_stencil,-1);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,0,o_stencil,1);
        
    }


;; {
;; _stencil_val(uf.x,1,o_stencil,0);_stencil_val(uf.x,0,o_stencil,0);         
       _stencil_val(f,-1,o_stencil,0);_stencil_val(f,0, o_stencil,0);_stencil_val(f,0, o_stencil,0); _stencil_val(f,1,o_stencil,0);
        
    }


_stencil_val(uf.y,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.y,0,0,0);  
  }}end__stencil_is_face_y()
#line 29
_stencil_is_face_z(){ {







;_stencil_val(uf.z,0,0,0);              
    
    _stencil_val(g.z,0,0,o_stencil);;;_stencil_val(f,0,0, o_stencil);





;; {
;; _stencil_val(uf.x,1,0,o_stencil);_stencil_val(uf.x,0,0,o_stencil);         
       _stencil_val(f,-1,0,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,1,0,o_stencil);
        
    }


;; {
;; _stencil_val(uf.y,0,1,o_stencil);_stencil_val(uf.y,0,0,o_stencil);         
       _stencil_val(f,0,-1,o_stencil);_stencil_val(f,0,0, o_stencil);_stencil_val(f,0,0, o_stencil); _stencil_val(f,0,1,o_stencil);
        
    }


_stencil_val(uf.z,0,0,0);





      


      


    
#line 59
_stencil_val_a(flux.z,0,0,0);  
  }}end__stencil_is_face_z()}end_foreach_face_stencil()




   BEGIN_FOREACH{
#line 29
foreach_face_generic(){is_face_x(){ {







    double un = dt*val(uf.x,0,0,0)/(_const_fm.x*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (_const_fm.y && _const_fm.y) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(_const_fm.y + _const_fm.y);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.z && _const_fm.z) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(_const_fm.z + _const_fm.z);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  }}end_is_face_x()
#line 29
is_face_y(){ {







    double un = dt*val(uf.y,0,0,0)/(_const_fm.y*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (_const_fm.z && _const_fm.z) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(_const_fm.z + _const_fm.z);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.x && _const_fm.x) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(_const_fm.x + _const_fm.x);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  }}end_is_face_y()
#line 29
is_face_z(){ {







    double un = dt*val(uf.z,0,0,0)/(_const_fm.z*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (_const_src + _const_src)*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (_const_fm.x && _const_fm.x) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(_const_fm.x + _const_fm.x);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (_const_fm.y && _const_fm.y) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(_const_fm.y + _const_fm.y);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  }}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }delete((scalar*)((vector[]){g,{{-1},{-1},{-1}}}));
end_tracing("tracer_fluxes","/home/esteban/ProgramFile/basilisk/src/bcg.h",61);}






     
void advection (scalar * tracers, vector u, double dt,
  scalar * src)
{tracing("advection","/home/esteban/ProgramFile/basilisk/src/bcg.h",69);




  scalar * psrc = src;
  if (!src)
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      const scalar zero = new_const_scalar("zero",16, 0.);
      src = list_append (src, zero);
    }}}
  if (!(list_len (tracers) == list_len (src))) qassert ("/home/esteban/ProgramFile/basilisk/src/bcg.h", 82, "list_len (tracers) == list_len (src)");

  scalar f, source;
  {scalar*_i0=src;scalar*_i1= tracers;if(_i0)for(source=*_i0,f=*_i1;_i0->i>= 0;source=*++_i0,f=*++_i1){ {
    vector  flux=new_face_vector("flux");
    tracer_fluxes (f, u, flux, dt, source);

    if(!is_constant(cm)){

    
#line 89
foreach_stencil(1,{0},)
      {
        {_stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);_stencil_val(cm,0,0,0);_stencil_val_r(f,0,0,0);   }
        
#line 91
{_stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);_stencil_val(cm,0,0,0);_stencil_val_r(f,0,0,0);   }
        
#line 91
{_stencil_val(flux.z,0,0,0); _stencil_val(flux.z,0,0,1);_stencil_val(cm,0,0,0);_stencil_val_r(f,0,0,0);   }}end_foreach_stencil() BEGIN_FOREACH{
#line 89
foreach()
      {
        val(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val(cm,0,0,0));
        
#line 91
val(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val(cm,0,0,0));
        
#line 91
val(f,0,0,0) += dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*val(cm,0,0,0));}end_foreach();}END_FOREACH }else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

    
#line 89
foreach_stencil(1,{0},)
      {
        {_stencil_val(flux.x,0,0,0); _stencil_val(flux.x,1,0,0);;_stencil_val_r(f,0,0,0);   }
        
#line 91
{_stencil_val(flux.y,0,0,0); _stencil_val(flux.y,0,1,0);;_stencil_val_r(f,0,0,0);   }
        
#line 91
{_stencil_val(flux.z,0,0,0); _stencil_val(flux.z,0,0,1);;_stencil_val_r(f,0,0,0);   }}end_foreach_stencil()

     BEGIN_FOREACH{
#line 89
foreach()
      {
        val(f,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*_const_cm);
        
#line 91
val(f,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*_const_cm);
        
#line 91
val(f,0,0,0) += dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*_const_cm);}end_foreach();}END_FOREACH }delete((scalar*)((vector[]){flux,{{-1},{-1},{-1}}}));



  }}}

  if (!psrc)
    pfree (src,__func__,__FILE__,__LINE__);
end_tracing("advection","/home/esteban/ProgramFile/basilisk/src/bcg.h",99);}
#line 30 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"



#line 1 "./viscosity.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
#line 50 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
#line 1 "./poisson.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
#line 32 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
     
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{tracing("mg_cycle","/home/esteban/ProgramFile/basilisk/src/poisson.h",33);




  restriction (res);





  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
      {
      
#line 56
foreach_level_or_leaf_stencil (1,DEPAREN({l}),)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     {_stencil_val_a(s,0,0,0);  }}}end_foreach_level_or_leaf_stencil() BEGIN_FOREACH{
#line 56
foreach_level_or_leaf (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = 0.;}}end_foreach_level_or_leaf();}END_FOREACH }





    else {
      boundary_level (da, l - 1);
      foreach_level_stencil (1,DEPAREN({l}),)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     { _stencil_bilinear (point, s);_stencil_val_a(s,0,0,0); }}}end_foreach_level_stencil()
       BEGIN_FOREACH{
#line 67
foreach_level (l)
 {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
  
     val(s,0,0,0) = bilinear (point, s);}}end_foreach_level();}END_FOREACH 
    }





    for (int i = 0; i < nrelax; i++) {
      boundary_level (da, l);
      relax (da, res, l, data);
    }
  }




  foreach_stencil(1,{0},) {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 { _stencil_val(ds,0,0,0);_stencil_val_r(s,0,0,0); }}}
  }end_foreach_stencil()




   BEGIN_FOREACH{
#line 86
foreach() {
    scalar s, ds;
    {scalar*_i0= da;scalar*_i1= a;if(_i0)for(ds=*_i0,s=*_i1;_i0->i>= 0;ds=*++_i0,s=*++_i1){
     
 val(s,0,0,0) += val(ds,0,0,0);}}
  }end_foreach();}END_FOREACH 
end_tracing("mg_cycle","/home/esteban/ProgramFile/basilisk/src/poisson.h",92);}
#line 104 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
  int minlevel;
} mgstats;
#line 127 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
     
mgstats mg_solve (scalar * a, scalar * b,
    double (* residual) (scalar * a, scalar * b, scalar * res,
           void * data),
    void (* relax) (scalar * da, scalar * res, int depth,
      void * data),
    void * data,
    int nrelax,
    scalar * res,
    int minlevel,
    double tolerance)
{tracing("mg_solve","/home/esteban/ProgramFile/basilisk/src/poisson.h",128);





  scalar * da = list_clone (a), * pres = res;
  if (!res)
    res = list_clone (b);






  for (int b = 0; b < nboundary; b++)
    {scalar*_i=(scalar*)( da);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];}}




  mgstats s = {0};
  double sum = 0.;
  scalar rhs = b[0];
  foreach_stencil (1,{0},)
    { _stencil_val(rhs,0,0,0); }end_foreach_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(+:sum)){
#line 163
foreach ()
    sum += val(rhs,0,0,0);end_foreach();mpi_all_reduce_array(&sum,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
  
#line 165
s.sum = sum;
  s.nrelax = nrelax > 0 ? nrelax : 4;




  double resb;
  resb = s.resb = s.resa = (* residual) (a, b, res, data);






  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > tolerance);
       s.i++) {
    mg_cycle (a, res, da, relax, data,
       s.nrelax,
       minlevel,
       grid->maxdepth);
    s.resa = (* residual) (a, b, res, data);
#line 195 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
    if (s.resa > tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }







    resb = s.resa;
  }
  s.minlevel = minlevel;




  if (s.resa > tolerance) {
    scalar v = a[0];
    fprintf (ferr,
      "src/poisson.h:%d: warning: convergence for %s not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d tolerance: %g\n", 219, _attribute[v.i].name,
      s.i, s.resa, s.sum, s.nrelax, tolerance), fflush (ferr);
  }




  if (!pres)
    delete (res), pfree (res,__func__,__FILE__,__LINE__);
  delete (da), pfree (da,__func__,__FILE__,__LINE__);

  {end_tracing("mg_solve","/home/esteban/ProgramFile/basilisk/src/poisson.h",230);return s;}
end_tracing("mg_solve","/home/esteban/ProgramFile/basilisk/src/poisson.h",231);}
#line 254 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
          vector alpha;
          scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;



};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
#line 292 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  scalar c = a;
#line 305 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  if(!is_constant(lambda) && !is_constant(alpha.x)){
#line 305 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  foreach_level_or_leaf_stencil (1,DEPAREN({l}),)

  {       





     _stencil_val(lambda,0,0,0);_stencil_val(b,0,0,0);
     { 
_stencil_val(alpha.x,1,0,0);_stencil_val(a,1,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(a,-1,0,0);
         _stencil_val(alpha.x,1,0,0); _stencil_val(alpha.x,0,0,0);
        
    
#line 317
} 
#line 314
{ 
_stencil_val(alpha.y,0,1,0);_stencil_val(a,0,1,0); _stencil_val(alpha.y,0,0,0);_stencil_val(a,0,-1,0);
         _stencil_val(alpha.y,0,1,0); _stencil_val(alpha.y,0,0,0);
        
    
#line 317
} 
#line 314
{ 
_stencil_val(alpha.z,0,0,1);_stencil_val(a,0,0,1); _stencil_val(alpha.z,0,0,0);_stencil_val(a,0,0,-1);
         _stencil_val(alpha.z,0,0,1); _stencil_val(alpha.z,0,0,0);
        
    
#line 317
}
#line 328 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
      _stencil_val_a(c,0,0,0);  
  }end_foreach_level_or_leaf_stencil() BEGIN_FOREACH{
#line 305
foreach_level_or_leaf (l)

  {





    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 314
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    } 
#line 314
{
      n += val(alpha.z,0,0,1)*val(a,0,0,1) + val(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val(alpha.z,0,0,1) + val(alpha.z,0,0,0);
    }
#line 328 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}END_FOREACH }else if(is_constant(lambda) && !is_constant(alpha.x)){double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);
#line 305 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  foreach_level_or_leaf_stencil (1,DEPAREN({l}),)

  {





;_stencil_val(b,0,0,0);
     { 
_stencil_val(alpha.x,1,0,0);_stencil_val(a,1,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(a,-1,0,0);
         _stencil_val(alpha.x,1,0,0); _stencil_val(alpha.x,0,0,0);
        
    
#line 317
} 
#line 314
{ 
_stencil_val(alpha.y,0,1,0);_stencil_val(a,0,1,0); _stencil_val(alpha.y,0,0,0);_stencil_val(a,0,-1,0);
         _stencil_val(alpha.y,0,1,0); _stencil_val(alpha.y,0,0,0);
        
    
#line 317
} 
#line 314
{ 
_stencil_val(alpha.z,0,0,1);_stencil_val(a,0,0,1); _stencil_val(alpha.z,0,0,0);_stencil_val(a,0,0,-1);
         _stencil_val(alpha.z,0,0,1); _stencil_val(alpha.z,0,0,0);
        
    
#line 317
}
#line 328 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
      _stencil_val_a(c,0,0,0);  
  }end_foreach_level_or_leaf_stencil()
#line 305 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
   BEGIN_FOREACH{foreach_level_or_leaf (l)

  {





    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += val(alpha.x,1,0,0)*val(a,1,0,0) + val(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val(alpha.x,1,0,0) + val(alpha.x,0,0,0);
    } 
#line 314
{
      n += val(alpha.y,0,1,0)*val(a,0,1,0) + val(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val(alpha.y,0,1,0) + val(alpha.y,0,0,0);
    } 
#line 314
{
      n += val(alpha.z,0,0,1)*val(a,0,0,1) + val(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val(alpha.z,0,0,1) + val(alpha.z,0,0,0);
    }
#line 328 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}END_FOREACH }else if(!is_constant(lambda) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
#line 305 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  foreach_level_or_leaf_stencil (1,DEPAREN({l}),)

  {       





     _stencil_val(lambda,0,0,0);_stencil_val(b,0,0,0);
     {
;_stencil_val(a,1,0,0);;_stencil_val(a,-1,0,0);
;;
        
    
#line 317
} 
#line 314
{
;_stencil_val(a,0,1,0);;_stencil_val(a,0,-1,0);
;;
        
    
#line 317
} 
#line 314
{
;_stencil_val(a,0,0,1);;_stencil_val(a,0,0,-1);
;;
        
    
#line 317
}
#line 328 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
      _stencil_val_a(c,0,0,0);  
  }end_foreach_level_or_leaf_stencil()
#line 305 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
   BEGIN_FOREACH{foreach_level_or_leaf (l)

  {





    double n = - sq(Delta)*val(b,0,0,0), d = - val(lambda,0,0,0)*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 314
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    } 
#line 314
{
      n += _const_alpha.z*val(a,0,0,1) + _const_alpha.z*val(a,0,0,-1);
      d += _const_alpha.z + _const_alpha.z;
    }
#line 328 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}END_FOREACH }else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
#line 305 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  foreach_level_or_leaf_stencil (1,DEPAREN({l}),)

  {





;_stencil_val(b,0,0,0);
     {
;_stencil_val(a,1,0,0);;_stencil_val(a,-1,0,0);
;;
        
    
#line 317
} 
#line 314
{
;_stencil_val(a,0,1,0);;_stencil_val(a,0,-1,0);
;;
        
    
#line 317
} 
#line 314
{
;_stencil_val(a,0,0,1);;_stencil_val(a,0,0,-1);
;;
        
    
#line 317
}
#line 328 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
      _stencil_val_a(c,0,0,0);  
  }end_foreach_level_or_leaf_stencil()
#line 305 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
   BEGIN_FOREACH{foreach_level_or_leaf (l)

  {





    double n = - sq(Delta)*val(b,0,0,0), d = - _const_lambda*sq(Delta);
     {
      n += _const_alpha.x*val(a,1,0,0) + _const_alpha.x*val(a,-1,0,0);
      d += _const_alpha.x + _const_alpha.x;
    } 
#line 314
{
      n += _const_alpha.y*val(a,0,1,0) + _const_alpha.y*val(a,0,-1,0);
      d += _const_alpha.y + _const_alpha.y;
    } 
#line 314
{
      n += _const_alpha.z*val(a,0,0,1) + _const_alpha.z*val(a,0,0,-1);
      d += _const_alpha.z + _const_alpha.z;
    }
#line 328 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  }end_foreach_level_or_leaf();}END_FOREACH }
#line 347 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
          vector alpha = p->alpha;
          scalar lambda = p->lambda;
  double maxres = 0.;
#line 381 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  if(!is_constant(lambda) && !is_constant(alpha.x)){
#line 381 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  foreach_stencil (1,{0},) { 
_stencil_val(b,0,0,0); _stencil_val(lambda,0,0,0);_stencil_val(a,0,0,0);
    
#line 382
_stencil_val_a(res,0,0,0);  
    
      {_stencil_val(alpha.x,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);
  _stencil_val(alpha.x,1,0,0);_stencil_val(a,1,0,0); _stencil_val(a,1 -1,0,0);
#line 384
_stencil_val_r(res,0,0,0);     
}
      
#line 384
{_stencil_val(alpha.y,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);
  _stencil_val(alpha.y,0,1,0);_stencil_val(a,0,1,0); _stencil_val(a,0,1 -1,0);
#line 384
_stencil_val_r(res,0,0,0);     
}
      
#line 384
{_stencil_val(alpha.z,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0,0,0 -1);
  _stencil_val(alpha.z,0,0,1);_stencil_val(a,0,0,1); _stencil_val(a,0,0,1 -1);
#line 384
_stencil_val_r(res,0,0,0);     
}






_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }






        
  
#line 394
}end_foreach_stencil()
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(max:maxres)){
#line 381
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);
    
      val(res,0,0,0) += (val(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 384
val(res,0,0,0) += (val(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;
      
#line 384
val(res,0,0,0) += (val(alpha.z,0,0,0)*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta) -
  val(alpha.z,0,0,1)*((val(a,0,0,1) - val(a,0,0,1 -1))/Delta))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 394
}else if(is_constant(lambda) && !is_constant(alpha.x)){double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);
#line 381 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  foreach_stencil (1,{0},) { 
_stencil_val(b,0,0,0);;_stencil_val(a,0,0,0);
    
#line 382
_stencil_val_a(res,0,0,0);  
    
      {_stencil_val(alpha.x,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);
  _stencil_val(alpha.x,1,0,0);_stencil_val(a,1,0,0); _stencil_val(a,1 -1,0,0);
#line 384
_stencil_val_r(res,0,0,0);     
}
      
#line 384
{_stencil_val(alpha.y,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);
  _stencil_val(alpha.y,0,1,0);_stencil_val(a,0,1,0); _stencil_val(a,0,1 -1,0);
#line 384
_stencil_val_r(res,0,0,0);     
}
      
#line 384
{_stencil_val(alpha.z,0,0,0);_stencil_val(a,0,0,0); _stencil_val(a,0,0,0 -1);
  _stencil_val(alpha.z,0,0,1);_stencil_val(a,0,0,1); _stencil_val(a,0,0,1 -1);
#line 384
_stencil_val_r(res,0,0,0);     
}






_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }






        
  
#line 394
}end_foreach_stencil()
#line 381 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(max:maxres)){
#line 381
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);
    
      val(res,0,0,0) += (val(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 384
val(res,0,0,0) += (val(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;
      
#line 384
val(res,0,0,0) += (val(alpha.z,0,0,0)*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta) -
  val(alpha.z,0,0,1)*((val(a,0,0,1) - val(a,0,0,1 -1))/Delta))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 394
}else if(!is_constant(lambda) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
#line 381 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  foreach_stencil (1,{0},) { 
_stencil_val(b,0,0,0); _stencil_val(lambda,0,0,0);_stencil_val(a,0,0,0);
    
#line 382
_stencil_val_a(res,0,0,0);  
    
      {;_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);
;_stencil_val(a,1,0,0); _stencil_val(a,1 -1,0,0);
#line 384
_stencil_val_r(res,0,0,0);     
}
      
#line 384
{;_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);
;_stencil_val(a,0,1,0); _stencil_val(a,0,1 -1,0);
#line 384
_stencil_val_r(res,0,0,0);     
}
      
#line 384
{;_stencil_val(a,0,0,0); _stencil_val(a,0,0,0 -1);
;_stencil_val(a,0,0,1); _stencil_val(a,0,0,1 -1);
#line 384
_stencil_val_r(res,0,0,0);     
}






_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }






        
  
#line 394
}end_foreach_stencil()
#line 381 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(max:maxres)){
#line 381
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - val(lambda,0,0,0)*val(a,0,0,0);
    
      val(res,0,0,0) += (_const_alpha.x*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  _const_alpha.x*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 384
val(res,0,0,0) += (_const_alpha.y*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  _const_alpha.y*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;
      
#line 384
val(res,0,0,0) += (_const_alpha.z*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta) -
  _const_alpha.z*((val(a,0,0,1) - val(a,0,0,1 -1))/Delta))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 394
}else {double _const_lambda=_constant[lambda.i-_NVARMAX];NOT_UNUSED(_const_lambda);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
#line 381 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  foreach_stencil (1,{0},) { 
_stencil_val(b,0,0,0);;_stencil_val(a,0,0,0);
    
#line 382
_stencil_val_a(res,0,0,0);  
    
      {;_stencil_val(a,0,0,0); _stencil_val(a,0 -1,0,0);
;_stencil_val(a,1,0,0); _stencil_val(a,1 -1,0,0);
#line 384
_stencil_val_r(res,0,0,0);     
}
      
#line 384
{;_stencil_val(a,0,0,0); _stencil_val(a,0,0 -1,0);
;_stencil_val(a,0,1,0); _stencil_val(a,0,1 -1,0);
#line 384
_stencil_val_r(res,0,0,0);     
}
      
#line 384
{;_stencil_val(a,0,0,0); _stencil_val(a,0,0,0 -1);
;_stencil_val(a,0,0,1); _stencil_val(a,0,0,1 -1);
#line 384
_stencil_val_r(res,0,0,0);     
}






_stencil_val(res,0,0,0);
      {_stencil_val(res,0,0,0);   }






        
  
#line 394
}end_foreach_stencil()
#line 381 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(max:maxres)){
#line 381
foreach () {
    val(res,0,0,0) = val(b,0,0,0) - _const_lambda*val(a,0,0,0);
    
      val(res,0,0,0) += (_const_alpha.x*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  _const_alpha.x*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
      
#line 384
val(res,0,0,0) += (_const_alpha.y*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  _const_alpha.y*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;
      
#line 384
val(res,0,0,0) += (_const_alpha.z*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta) -
  _const_alpha.z*((val(a,0,0,1) - val(a,0,0,1 -1))/Delta))/Delta;






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  }end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 394
}

  return maxres;
}
#line 408 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
     
mgstats poisson (scalar a, scalar b,
           vector alpha,
           scalar lambda,
   double tolerance,
   int nrelax,
   int minlevel,
   scalar * res,
   double (* flux) (Point, scalar, vector, double *))
{tracing("poisson","/home/esteban/ProgramFile/basilisk/src/poisson.h",409);






  if (alpha.x.i < 0)
    alpha = new_const_vector("alpha",0,(double[]) {1.,1.,1.});
  if (lambda.i < 0)
    lambda = new_const_scalar("lambda",3, 0.);






  restriction (((scalar[]){alpha.x,alpha.y,alpha.z,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (tolerance)
    TOLERANCE = tolerance;

  struct Poisson p = {a, b, alpha, lambda, tolerance, nrelax, minlevel, res };






  mgstats s = mg_solve ((
#line 128
scalar *
#line 451
)((scalar[]){a,{-1}}),( 
#line 128
scalar *
#line 451
)((scalar[]) {b,{-1}}), residual, relax, &p
,
   
#line 452
nrelax, res, max(1, minlevel)
#line 136
, 
TOLERANCE
#line 452
);




  if (tolerance)
    TOLERANCE = defaultol;

  {end_tracing("poisson","/home/esteban/ProgramFile/basilisk/src/poisson.h",460);return s;}
end_tracing("poisson","/home/esteban/ProgramFile/basilisk/src/poisson.h",461);}
#line 480 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
     
mgstats project (vector uf, scalar p,
           vector alpha,
   double dt,
   int nrelax)
{tracing("project","/home/esteban/ProgramFile/basilisk/src/poisson.h",481);






  scalar  div=new_scalar("div");
  foreach_stencil(1,{0},) {
    _stencil_val_a(div,0,0,0);  
    
      { _stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);_stencil_val_r(div,0,0,0);  }
      
#line 496
{ _stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);_stencil_val_r(div,0,0,0);  }
      
#line 496
{ _stencil_val(uf.z,0,0,1); _stencil_val(uf.z,0,0,0);_stencil_val_r(div,0,0,0);  }
    _stencil_val_r(div,0,0,0);  
  }end_foreach_stencil()
   BEGIN_FOREACH{
#line 493
foreach() {
    val(div,0,0,0) = 0.;
    
      val(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);
      
#line 496
val(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);
      
#line 496
val(div,0,0,0) += val(uf.z,0,0,1) - val(uf.z,0,0,0);
    val(div,0,0,0) /= dt*Delta;
  }end_foreach();}END_FOREACH 
#line 509 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
  mgstats mgp = poisson (p, div, alpha
#line 410
,
( scalar) {-1}
#line 510
, TOLERANCE/sq(dt), nrelax
#line 413
, 
0, 
NULL, 
NULL
#line 510
);




  if(!is_constant(alpha.x)){




  
#line 515
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {_stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0 -1,0,0);_stencil_val_r(uf.x,0,0,0);   }}end__stencil_is_face_x()
#line 515
_stencil_is_face_y(){
    {_stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0 -1,0);_stencil_val_r(uf.y,0,0,0);   }}end__stencil_is_face_y()
#line 515
_stencil_is_face_z(){
    {_stencil_val(alpha.z,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0,0 -1);_stencil_val_r(uf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil() BEGIN_FOREACH{
#line 515
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*val(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 515
is_face_y(){
    val(uf.y,0,0,0) -= dt*val(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()
#line 515
is_face_z(){
    val(uf.z,0,0,0) -= dt*val(alpha.z,0,0,0)*((val(p,0,0,0) - val(p,0,0,0 -1))/Delta);}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else {_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);




  
#line 515
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {;_stencil_val(p,0,0,0); _stencil_val(p,0 -1,0,0);_stencil_val_r(uf.x,0,0,0);   }}end__stencil_is_face_x()
#line 515
_stencil_is_face_y(){
    {;_stencil_val(p,0,0,0); _stencil_val(p,0,0 -1,0);_stencil_val_r(uf.y,0,0,0);   }}end__stencil_is_face_y()
#line 515
_stencil_is_face_z(){
    {;_stencil_val(p,0,0,0); _stencil_val(p,0,0,0 -1);_stencil_val_r(uf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil()




   BEGIN_FOREACH{
#line 515
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) -= dt*_const_alpha.x*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta);}end_is_face_x()
#line 515
is_face_y(){
    val(uf.y,0,0,0) -= dt*_const_alpha.y*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta);}end_is_face_y()
#line 515
is_face_z(){
    val(uf.z,0,0,0) -= dt*_const_alpha.z*((val(p,0,0,0) - val(p,0,0,0 -1))/Delta);}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }

  {delete((scalar*)((scalar[]){div,{-1}}));{end_tracing("project","/home/esteban/ProgramFile/basilisk/src/poisson.h",518);return mgp;}}delete((scalar*)((scalar[]){div,{-1}}));
end_tracing("project","/home/esteban/ProgramFile/basilisk/src/poisson.h",519);}
#line 51 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"

struct Viscosity {
  vector mu;
  scalar rho;
  double dt;
};
#line 129 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));




  vector w = u;
#line 161 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
  vector ua = u;

  if(!is_constant(rho) && !is_constant(mu.x)){

  
#line 163
foreach_level_or_leaf_stencil (1,DEPAREN({l}),)

  {
    
      {_stencil_val(r.x,0,0,0);_stencil_val(rho,0,0,0);_stencil_val(mu.x,1,0,0);_stencil_val(u.x,1,0,0);_stencil_val(mu.x,0,0,0);_stencil_val(u.x,-1,0,0); 

_stencil_val(mu.y,0,1,0);_stencil_val(u.x,0,1,0);
_stencil_val(u.y,1,0,0); _stencil_val(ua.y,1,1,0);
_stencil_val(u.y,-1,0,0); _stencil_val(ua.y,-1,1,0); 
_stencil_val(mu.y,0,0,0); _stencil_val(u.x,0,-1,0);
_stencil_val(ua.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(ua.y,-1,-1,0); _stencil_val(u.y,-1,0,0); 


_stencil_val(mu.z,0,0,1);_stencil_val(u.x,0,0,1);
_stencil_val(u.z,1,0,0); _stencil_val(ua.z,1,0,1);
_stencil_val(u.z,-1,0,0); _stencil_val(ua.z,-1,0,1); 
_stencil_val(mu.z,0,0,0); _stencil_val(u.x,0,0,-1);
_stencil_val(ua.z,1,0,-1); _stencil_val(u.z,1,0,0);
_stencil_val(ua.z,-1,0,-1); _stencil_val(u.z,-1,0,0);


_stencil_val(rho,0,0,0);_stencil_val(mu.x,1,0,0);_stencil_val(mu.x,0,0,0); 

_stencil_val(mu.y,0,1,0); _stencil_val(mu.y,0,0,0); 


_stencil_val(mu.z,0,0,1); _stencil_val(mu.z,0,0,0);
#line 167
_stencil_val_a(w.x,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
      
#line 167
{_stencil_val(r.y,0,0,0);_stencil_val(rho,0,0,0);_stencil_val(mu.y,0,1,0);_stencil_val(u.y,0,1,0);_stencil_val(mu.y,0,0,0);_stencil_val(u.y,0,-1,0); 

_stencil_val(mu.z,0,0,1);_stencil_val(u.y,0,0,1);
_stencil_val(u.z,0,1,0); _stencil_val(ua.z,0,1,1);
_stencil_val(u.z,0,-1,0); _stencil_val(ua.z,0,-1,1); 
_stencil_val(mu.z,0,0,0); _stencil_val(u.y,0,0,-1);
_stencil_val(ua.z,0,1,-1); _stencil_val(u.z,0,1,0);
_stencil_val(ua.z,0,-1,-1); _stencil_val(u.z,0,-1,0); 


_stencil_val(mu.x,1,0,0);_stencil_val(u.y,1,0,0);
_stencil_val(u.x,0,1,0); _stencil_val(ua.x,1,1,0);
_stencil_val(u.x,0,-1,0); _stencil_val(ua.x,1,-1,0); 
_stencil_val(mu.x,0,0,0); _stencil_val(u.y,-1,0,0);
_stencil_val(ua.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(ua.x,-1,-1,0); _stencil_val(u.x,0,-1,0);


_stencil_val(rho,0,0,0);_stencil_val(mu.y,0,1,0);_stencil_val(mu.y,0,0,0); 

_stencil_val(mu.z,0,0,1); _stencil_val(mu.z,0,0,0); 


_stencil_val(mu.x,1,0,0); _stencil_val(mu.x,0,0,0);
#line 167
_stencil_val_a(w.y,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
      
#line 167
{_stencil_val(r.z,0,0,0);_stencil_val(rho,0,0,0);_stencil_val(mu.z,0,0,1);_stencil_val(u.z,0,0,1);_stencil_val(mu.z,0,0,0);_stencil_val(u.z,0,0,-1); 

_stencil_val(mu.x,1,0,0);_stencil_val(u.z,1,0,0);
_stencil_val(u.x,0,0,1); _stencil_val(ua.x,1,0,1);
_stencil_val(u.x,0,0,-1); _stencil_val(ua.x,1,0,-1); 
_stencil_val(mu.x,0,0,0); _stencil_val(u.z,-1,0,0);
_stencil_val(ua.x,-1,0,1); _stencil_val(u.x,0,0,1);
_stencil_val(ua.x,-1,0,-1); _stencil_val(u.x,0,0,-1); 


_stencil_val(mu.y,0,1,0);_stencil_val(u.z,0,1,0);
_stencil_val(u.y,0,0,1); _stencil_val(ua.y,0,1,1);
_stencil_val(u.y,0,0,-1); _stencil_val(ua.y,0,1,-1); 
_stencil_val(mu.y,0,0,0); _stencil_val(u.z,0,-1,0);
_stencil_val(ua.y,0,-1,1); _stencil_val(u.y,0,0,1);
_stencil_val(ua.y,0,-1,-1); _stencil_val(u.y,0,0,-1);


_stencil_val(rho,0,0,0);_stencil_val(mu.z,0,0,1);_stencil_val(mu.z,0,0,0); 

_stencil_val(mu.x,1,0,0); _stencil_val(mu.x,0,0,0); 


_stencil_val(mu.y,0,1,0); _stencil_val(mu.y,0,0,0);
#line 167
_stencil_val_a(w.z,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
  }end_foreach_level_or_leaf_stencil() BEGIN_FOREACH{
#line 163
foreach_level_or_leaf (l)

  {
    
      val(w.x,0,0,0) = (val(r.x,0,0,0)*sq(Delta) + dt/val(rho,0,0,0)*(2.*val(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val(mu.x,0,0,0)*val(u.x,-1,0,0)

        + val(mu.y,0,1,0)*(val(u.x,0,1,0) +
       (val(u.y,1,0,0) + val(ua.y,1,1,0))/4. -
       (val(u.y,-1,0,0) + val(ua.y,-1,1,0))/4.)
        - val(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
           (val(ua.y,1,-1,0) + val(u.y,1,0,0))/4. -
           (val(ua.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


        + val(mu.z,0,0,1)*(val(u.x,0,0,1) +
         (val(u.z,1,0,0) + val(ua.z,1,0,1))/4. -
         (val(u.z,-1,0,0) + val(ua.z,-1,0,1))/4.)
        - val(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
           (val(ua.z,1,0,-1) + val(u.z,1,0,0))/4. -
           (val(ua.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

        ))/
      (((coord){1.,1.,1.}).x*sq(Delta) + dt/val(rho,0,0,0)*(2.*val(mu.x,1,0,0) + 2.*val(mu.x,0,0,0)

          + val(mu.y,0,1,0) + val(mu.y,0,0,0)


          + val(mu.z,0,0,1) + val(mu.z,0,0,0)

          ));
      
#line 167
val(w.y,0,0,0) = (val(r.y,0,0,0)*sq(Delta) + dt/val(rho,0,0,0)*(2.*val(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val(mu.y,0,0,0)*val(u.y,0,-1,0)

        + val(mu.z,0,0,1)*(val(u.y,0,0,1) +
       (val(u.z,0,1,0) + val(ua.z,0,1,1))/4. -
       (val(u.z,0,-1,0) + val(ua.z,0,-1,1))/4.)
        - val(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
           (val(ua.z,0,1,-1) + val(u.z,0,1,0))/4. -
           (val(ua.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


        + val(mu.x,1,0,0)*(val(u.y,1,0,0) +
         (val(u.x,0,1,0) + val(ua.x,1,1,0))/4. -
         (val(u.x,0,-1,0) + val(ua.x,1,-1,0))/4.)
        - val(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
           (val(ua.x,-1,1,0) + val(u.x,0,1,0))/4. -
           (val(ua.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

        ))/
      (((coord){1.,1.,1.}).y*sq(Delta) + dt/val(rho,0,0,0)*(2.*val(mu.y,0,1,0) + 2.*val(mu.y,0,0,0)

          + val(mu.z,0,0,1) + val(mu.z,0,0,0)


          + val(mu.x,1,0,0) + val(mu.x,0,0,0)

          ));
      
#line 167
val(w.z,0,0,0) = (val(r.z,0,0,0)*sq(Delta) + dt/val(rho,0,0,0)*(2.*val(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val(mu.z,0,0,0)*val(u.z,0,0,-1)

        + val(mu.x,1,0,0)*(val(u.z,1,0,0) +
       (val(u.x,0,0,1) + val(ua.x,1,0,1))/4. -
       (val(u.x,0,0,-1) + val(ua.x,1,0,-1))/4.)
        - val(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
           (val(ua.x,-1,0,1) + val(u.x,0,0,1))/4. -
           (val(ua.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


        + val(mu.y,0,1,0)*(val(u.z,0,1,0) +
         (val(u.y,0,0,1) + val(ua.y,0,1,1))/4. -
         (val(u.y,0,0,-1) + val(ua.y,0,1,-1))/4.)
        - val(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
           (val(ua.y,0,-1,1) + val(u.y,0,0,1))/4. -
           (val(ua.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

        ))/
      (((coord){1.,1.,1.}).z*sq(Delta) + dt/val(rho,0,0,0)*(2.*val(mu.z,0,0,1) + 2.*val(mu.z,0,0,0)

          + val(mu.x,1,0,0) + val(mu.x,0,0,0)


          + val(mu.y,0,1,0) + val(mu.y,0,0,0)

          ));
  }end_foreach_level_or_leaf();}END_FOREACH }else if(is_constant(rho) && !is_constant(mu.x)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);

  
#line 163
foreach_level_or_leaf_stencil (1,DEPAREN({l}),)

  {
    
      {_stencil_val(r.x,0,0,0);;_stencil_val(mu.x,1,0,0);_stencil_val(u.x,1,0,0);_stencil_val(mu.x,0,0,0);_stencil_val(u.x,-1,0,0); 

_stencil_val(mu.y,0,1,0);_stencil_val(u.x,0,1,0);
_stencil_val(u.y,1,0,0); _stencil_val(ua.y,1,1,0);
_stencil_val(u.y,-1,0,0); _stencil_val(ua.y,-1,1,0); 
_stencil_val(mu.y,0,0,0); _stencil_val(u.x,0,-1,0);
_stencil_val(ua.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(ua.y,-1,-1,0); _stencil_val(u.y,-1,0,0); 


_stencil_val(mu.z,0,0,1);_stencil_val(u.x,0,0,1);
_stencil_val(u.z,1,0,0); _stencil_val(ua.z,1,0,1);
_stencil_val(u.z,-1,0,0); _stencil_val(ua.z,-1,0,1); 
_stencil_val(mu.z,0,0,0); _stencil_val(u.x,0,0,-1);
_stencil_val(ua.z,1,0,-1); _stencil_val(u.z,1,0,0);
_stencil_val(ua.z,-1,0,-1); _stencil_val(u.z,-1,0,0);


;_stencil_val(mu.x,1,0,0);_stencil_val(mu.x,0,0,0); 

_stencil_val(mu.y,0,1,0); _stencil_val(mu.y,0,0,0); 


_stencil_val(mu.z,0,0,1); _stencil_val(mu.z,0,0,0);
#line 167
_stencil_val_a(w.x,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
      
#line 167
{_stencil_val(r.y,0,0,0);;_stencil_val(mu.y,0,1,0);_stencil_val(u.y,0,1,0);_stencil_val(mu.y,0,0,0);_stencil_val(u.y,0,-1,0); 

_stencil_val(mu.z,0,0,1);_stencil_val(u.y,0,0,1);
_stencil_val(u.z,0,1,0); _stencil_val(ua.z,0,1,1);
_stencil_val(u.z,0,-1,0); _stencil_val(ua.z,0,-1,1); 
_stencil_val(mu.z,0,0,0); _stencil_val(u.y,0,0,-1);
_stencil_val(ua.z,0,1,-1); _stencil_val(u.z,0,1,0);
_stencil_val(ua.z,0,-1,-1); _stencil_val(u.z,0,-1,0); 


_stencil_val(mu.x,1,0,0);_stencil_val(u.y,1,0,0);
_stencil_val(u.x,0,1,0); _stencil_val(ua.x,1,1,0);
_stencil_val(u.x,0,-1,0); _stencil_val(ua.x,1,-1,0); 
_stencil_val(mu.x,0,0,0); _stencil_val(u.y,-1,0,0);
_stencil_val(ua.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(ua.x,-1,-1,0); _stencil_val(u.x,0,-1,0);


;_stencil_val(mu.y,0,1,0);_stencil_val(mu.y,0,0,0); 

_stencil_val(mu.z,0,0,1); _stencil_val(mu.z,0,0,0); 


_stencil_val(mu.x,1,0,0); _stencil_val(mu.x,0,0,0);
#line 167
_stencil_val_a(w.y,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
      
#line 167
{_stencil_val(r.z,0,0,0);;_stencil_val(mu.z,0,0,1);_stencil_val(u.z,0,0,1);_stencil_val(mu.z,0,0,0);_stencil_val(u.z,0,0,-1); 

_stencil_val(mu.x,1,0,0);_stencil_val(u.z,1,0,0);
_stencil_val(u.x,0,0,1); _stencil_val(ua.x,1,0,1);
_stencil_val(u.x,0,0,-1); _stencil_val(ua.x,1,0,-1); 
_stencil_val(mu.x,0,0,0); _stencil_val(u.z,-1,0,0);
_stencil_val(ua.x,-1,0,1); _stencil_val(u.x,0,0,1);
_stencil_val(ua.x,-1,0,-1); _stencil_val(u.x,0,0,-1); 


_stencil_val(mu.y,0,1,0);_stencil_val(u.z,0,1,0);
_stencil_val(u.y,0,0,1); _stencil_val(ua.y,0,1,1);
_stencil_val(u.y,0,0,-1); _stencil_val(ua.y,0,1,-1); 
_stencil_val(mu.y,0,0,0); _stencil_val(u.z,0,-1,0);
_stencil_val(ua.y,0,-1,1); _stencil_val(u.y,0,0,1);
_stencil_val(ua.y,0,-1,-1); _stencil_val(u.y,0,0,-1);


;_stencil_val(mu.z,0,0,1);_stencil_val(mu.z,0,0,0); 

_stencil_val(mu.x,1,0,0); _stencil_val(mu.x,0,0,0); 


_stencil_val(mu.y,0,1,0); _stencil_val(mu.y,0,0,0);
#line 167
_stencil_val_a(w.z,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
  }end_foreach_level_or_leaf_stencil()

   BEGIN_FOREACH{
#line 163
foreach_level_or_leaf (l)

  {
    
      val(w.x,0,0,0) = (val(r.x,0,0,0)*sq(Delta) + dt/_const_rho*(2.*val(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val(mu.x,0,0,0)*val(u.x,-1,0,0)

        + val(mu.y,0,1,0)*(val(u.x,0,1,0) +
       (val(u.y,1,0,0) + val(ua.y,1,1,0))/4. -
       (val(u.y,-1,0,0) + val(ua.y,-1,1,0))/4.)
        - val(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
           (val(ua.y,1,-1,0) + val(u.y,1,0,0))/4. -
           (val(ua.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


        + val(mu.z,0,0,1)*(val(u.x,0,0,1) +
         (val(u.z,1,0,0) + val(ua.z,1,0,1))/4. -
         (val(u.z,-1,0,0) + val(ua.z,-1,0,1))/4.)
        - val(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
           (val(ua.z,1,0,-1) + val(u.z,1,0,0))/4. -
           (val(ua.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

        ))/
      (((coord){1.,1.,1.}).x*sq(Delta) + dt/_const_rho*(2.*val(mu.x,1,0,0) + 2.*val(mu.x,0,0,0)

          + val(mu.y,0,1,0) + val(mu.y,0,0,0)


          + val(mu.z,0,0,1) + val(mu.z,0,0,0)

          ));
      
#line 167
val(w.y,0,0,0) = (val(r.y,0,0,0)*sq(Delta) + dt/_const_rho*(2.*val(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val(mu.y,0,0,0)*val(u.y,0,-1,0)

        + val(mu.z,0,0,1)*(val(u.y,0,0,1) +
       (val(u.z,0,1,0) + val(ua.z,0,1,1))/4. -
       (val(u.z,0,-1,0) + val(ua.z,0,-1,1))/4.)
        - val(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
           (val(ua.z,0,1,-1) + val(u.z,0,1,0))/4. -
           (val(ua.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


        + val(mu.x,1,0,0)*(val(u.y,1,0,0) +
         (val(u.x,0,1,0) + val(ua.x,1,1,0))/4. -
         (val(u.x,0,-1,0) + val(ua.x,1,-1,0))/4.)
        - val(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
           (val(ua.x,-1,1,0) + val(u.x,0,1,0))/4. -
           (val(ua.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

        ))/
      (((coord){1.,1.,1.}).y*sq(Delta) + dt/_const_rho*(2.*val(mu.y,0,1,0) + 2.*val(mu.y,0,0,0)

          + val(mu.z,0,0,1) + val(mu.z,0,0,0)


          + val(mu.x,1,0,0) + val(mu.x,0,0,0)

          ));
      
#line 167
val(w.z,0,0,0) = (val(r.z,0,0,0)*sq(Delta) + dt/_const_rho*(2.*val(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val(mu.z,0,0,0)*val(u.z,0,0,-1)

        + val(mu.x,1,0,0)*(val(u.z,1,0,0) +
       (val(u.x,0,0,1) + val(ua.x,1,0,1))/4. -
       (val(u.x,0,0,-1) + val(ua.x,1,0,-1))/4.)
        - val(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
           (val(ua.x,-1,0,1) + val(u.x,0,0,1))/4. -
           (val(ua.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


        + val(mu.y,0,1,0)*(val(u.z,0,1,0) +
         (val(u.y,0,0,1) + val(ua.y,0,1,1))/4. -
         (val(u.y,0,0,-1) + val(ua.y,0,1,-1))/4.)
        - val(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
           (val(ua.y,0,-1,1) + val(u.y,0,0,1))/4. -
           (val(ua.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

        ))/
      (((coord){1.,1.,1.}).z*sq(Delta) + dt/_const_rho*(2.*val(mu.z,0,0,1) + 2.*val(mu.z,0,0,0)

          + val(mu.x,1,0,0) + val(mu.x,0,0,0)


          + val(mu.y,0,1,0) + val(mu.y,0,0,0)

          ));
  }end_foreach_level_or_leaf();}END_FOREACH }else if(!is_constant(rho) && is_constant(mu.x)){_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX]};NOT_UNUSED(_const_mu);

  
#line 163
foreach_level_or_leaf_stencil (1,DEPAREN({l}),)

  {
    
      {_stencil_val(r.x,0,0,0);_stencil_val(rho,0,0,0);;_stencil_val(u.x,1,0,0);;_stencil_val(u.x,-1,0,0);

;_stencil_val(u.x,0,1,0);
_stencil_val(u.y,1,0,0); _stencil_val(ua.y,1,1,0);
_stencil_val(u.y,-1,0,0); _stencil_val(ua.y,-1,1,0);
; _stencil_val(u.x,0,-1,0);
_stencil_val(ua.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(ua.y,-1,-1,0); _stencil_val(u.y,-1,0,0);


;_stencil_val(u.x,0,0,1);
_stencil_val(u.z,1,0,0); _stencil_val(ua.z,1,0,1);
_stencil_val(u.z,-1,0,0); _stencil_val(ua.z,-1,0,1);
; _stencil_val(u.x,0,0,-1);
_stencil_val(ua.z,1,0,-1); _stencil_val(u.z,1,0,0);
_stencil_val(ua.z,-1,0,-1); _stencil_val(u.z,-1,0,0);


_stencil_val(rho,0,0,0);;;

;;


;;
#line 167
_stencil_val_a(w.x,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
      
#line 167
{_stencil_val(r.y,0,0,0);_stencil_val(rho,0,0,0);;_stencil_val(u.y,0,1,0);;_stencil_val(u.y,0,-1,0);

;_stencil_val(u.y,0,0,1);
_stencil_val(u.z,0,1,0); _stencil_val(ua.z,0,1,1);
_stencil_val(u.z,0,-1,0); _stencil_val(ua.z,0,-1,1);
; _stencil_val(u.y,0,0,-1);
_stencil_val(ua.z,0,1,-1); _stencil_val(u.z,0,1,0);
_stencil_val(ua.z,0,-1,-1); _stencil_val(u.z,0,-1,0);


;_stencil_val(u.y,1,0,0);
_stencil_val(u.x,0,1,0); _stencil_val(ua.x,1,1,0);
_stencil_val(u.x,0,-1,0); _stencil_val(ua.x,1,-1,0);
; _stencil_val(u.y,-1,0,0);
_stencil_val(ua.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(ua.x,-1,-1,0); _stencil_val(u.x,0,-1,0);


_stencil_val(rho,0,0,0);;;

;;


;;
#line 167
_stencil_val_a(w.y,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
      
#line 167
{_stencil_val(r.z,0,0,0);_stencil_val(rho,0,0,0);;_stencil_val(u.z,0,0,1);;_stencil_val(u.z,0,0,-1);

;_stencil_val(u.z,1,0,0);
_stencil_val(u.x,0,0,1); _stencil_val(ua.x,1,0,1);
_stencil_val(u.x,0,0,-1); _stencil_val(ua.x,1,0,-1);
; _stencil_val(u.z,-1,0,0);
_stencil_val(ua.x,-1,0,1); _stencil_val(u.x,0,0,1);
_stencil_val(ua.x,-1,0,-1); _stencil_val(u.x,0,0,-1);


;_stencil_val(u.z,0,1,0);
_stencil_val(u.y,0,0,1); _stencil_val(ua.y,0,1,1);
_stencil_val(u.y,0,0,-1); _stencil_val(ua.y,0,1,-1);
; _stencil_val(u.z,0,-1,0);
_stencil_val(ua.y,0,-1,1); _stencil_val(u.y,0,0,1);
_stencil_val(ua.y,0,-1,-1); _stencil_val(u.y,0,0,-1);


_stencil_val(rho,0,0,0);;;

;;


;;
#line 167
_stencil_val_a(w.z,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
  }end_foreach_level_or_leaf_stencil()

   BEGIN_FOREACH{
#line 163
foreach_level_or_leaf (l)

  {
    
      val(w.x,0,0,0) = (val(r.x,0,0,0)*sq(Delta) + dt/val(rho,0,0,0)*(2.*_const_mu.x*val(u.x,1,0,0) + 2.*_const_mu.x*val(u.x,-1,0,0)

        + _const_mu.y*(val(u.x,0,1,0) +
       (val(u.y,1,0,0) + val(ua.y,1,1,0))/4. -
       (val(u.y,-1,0,0) + val(ua.y,-1,1,0))/4.)
        - _const_mu.y*(- val(u.x,0,-1,0) +
           (val(ua.y,1,-1,0) + val(u.y,1,0,0))/4. -
           (val(ua.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


        + _const_mu.z*(val(u.x,0,0,1) +
         (val(u.z,1,0,0) + val(ua.z,1,0,1))/4. -
         (val(u.z,-1,0,0) + val(ua.z,-1,0,1))/4.)
        - _const_mu.z*(- val(u.x,0,0,-1) +
           (val(ua.z,1,0,-1) + val(u.z,1,0,0))/4. -
           (val(ua.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

        ))/
      (((coord){1.,1.,1.}).x*sq(Delta) + dt/val(rho,0,0,0)*(2.*_const_mu.x + 2.*_const_mu.x

          + _const_mu.y + _const_mu.y


          + _const_mu.z + _const_mu.z

          ));
      
#line 167
val(w.y,0,0,0) = (val(r.y,0,0,0)*sq(Delta) + dt/val(rho,0,0,0)*(2.*_const_mu.y*val(u.y,0,1,0) + 2.*_const_mu.y*val(u.y,0,-1,0)

        + _const_mu.z*(val(u.y,0,0,1) +
       (val(u.z,0,1,0) + val(ua.z,0,1,1))/4. -
       (val(u.z,0,-1,0) + val(ua.z,0,-1,1))/4.)
        - _const_mu.z*(- val(u.y,0,0,-1) +
           (val(ua.z,0,1,-1) + val(u.z,0,1,0))/4. -
           (val(ua.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


        + _const_mu.x*(val(u.y,1,0,0) +
         (val(u.x,0,1,0) + val(ua.x,1,1,0))/4. -
         (val(u.x,0,-1,0) + val(ua.x,1,-1,0))/4.)
        - _const_mu.x*(- val(u.y,-1,0,0) +
           (val(ua.x,-1,1,0) + val(u.x,0,1,0))/4. -
           (val(ua.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

        ))/
      (((coord){1.,1.,1.}).y*sq(Delta) + dt/val(rho,0,0,0)*(2.*_const_mu.y + 2.*_const_mu.y

          + _const_mu.z + _const_mu.z


          + _const_mu.x + _const_mu.x

          ));
      
#line 167
val(w.z,0,0,0) = (val(r.z,0,0,0)*sq(Delta) + dt/val(rho,0,0,0)*(2.*_const_mu.z*val(u.z,0,0,1) + 2.*_const_mu.z*val(u.z,0,0,-1)

        + _const_mu.x*(val(u.z,1,0,0) +
       (val(u.x,0,0,1) + val(ua.x,1,0,1))/4. -
       (val(u.x,0,0,-1) + val(ua.x,1,0,-1))/4.)
        - _const_mu.x*(- val(u.z,-1,0,0) +
           (val(ua.x,-1,0,1) + val(u.x,0,0,1))/4. -
           (val(ua.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


        + _const_mu.y*(val(u.z,0,1,0) +
         (val(u.y,0,0,1) + val(ua.y,0,1,1))/4. -
         (val(u.y,0,0,-1) + val(ua.y,0,1,-1))/4.)
        - _const_mu.y*(- val(u.z,0,-1,0) +
           (val(ua.y,0,-1,1) + val(u.y,0,0,1))/4. -
           (val(ua.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

        ))/
      (((coord){1.,1.,1.}).z*sq(Delta) + dt/val(rho,0,0,0)*(2.*_const_mu.z + 2.*_const_mu.z

          + _const_mu.x + _const_mu.x


          + _const_mu.y + _const_mu.y

          ));
  }end_foreach_level_or_leaf();}END_FOREACH }else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX]};NOT_UNUSED(_const_mu);

  
#line 163
foreach_level_or_leaf_stencil (1,DEPAREN({l}),)

  {
    
      {_stencil_val(r.x,0,0,0);;;_stencil_val(u.x,1,0,0);;_stencil_val(u.x,-1,0,0);

;_stencil_val(u.x,0,1,0);
_stencil_val(u.y,1,0,0); _stencil_val(ua.y,1,1,0);
_stencil_val(u.y,-1,0,0); _stencil_val(ua.y,-1,1,0);
; _stencil_val(u.x,0,-1,0);
_stencil_val(ua.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(ua.y,-1,-1,0); _stencil_val(u.y,-1,0,0);


;_stencil_val(u.x,0,0,1);
_stencil_val(u.z,1,0,0); _stencil_val(ua.z,1,0,1);
_stencil_val(u.z,-1,0,0); _stencil_val(ua.z,-1,0,1);
; _stencil_val(u.x,0,0,-1);
_stencil_val(ua.z,1,0,-1); _stencil_val(u.z,1,0,0);
_stencil_val(ua.z,-1,0,-1); _stencil_val(u.z,-1,0,0);


;;;

;;


;;
#line 167
_stencil_val_a(w.x,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
      
#line 167
{_stencil_val(r.y,0,0,0);;;_stencil_val(u.y,0,1,0);;_stencil_val(u.y,0,-1,0);

;_stencil_val(u.y,0,0,1);
_stencil_val(u.z,0,1,0); _stencil_val(ua.z,0,1,1);
_stencil_val(u.z,0,-1,0); _stencil_val(ua.z,0,-1,1);
; _stencil_val(u.y,0,0,-1);
_stencil_val(ua.z,0,1,-1); _stencil_val(u.z,0,1,0);
_stencil_val(ua.z,0,-1,-1); _stencil_val(u.z,0,-1,0);


;_stencil_val(u.y,1,0,0);
_stencil_val(u.x,0,1,0); _stencil_val(ua.x,1,1,0);
_stencil_val(u.x,0,-1,0); _stencil_val(ua.x,1,-1,0);
; _stencil_val(u.y,-1,0,0);
_stencil_val(ua.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(ua.x,-1,-1,0); _stencil_val(u.x,0,-1,0);


;;;

;;


;;
#line 167
_stencil_val_a(w.y,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
      
#line 167
{_stencil_val(r.z,0,0,0);;;_stencil_val(u.z,0,0,1);;_stencil_val(u.z,0,0,-1);

;_stencil_val(u.z,1,0,0);
_stencil_val(u.x,0,0,1); _stencil_val(ua.x,1,0,1);
_stencil_val(u.x,0,0,-1); _stencil_val(ua.x,1,0,-1);
; _stencil_val(u.z,-1,0,0);
_stencil_val(ua.x,-1,0,1); _stencil_val(u.x,0,0,1);
_stencil_val(ua.x,-1,0,-1); _stencil_val(u.x,0,0,-1);


;_stencil_val(u.z,0,1,0);
_stencil_val(u.y,0,0,1); _stencil_val(ua.y,0,1,1);
_stencil_val(u.y,0,0,-1); _stencil_val(ua.y,0,1,-1);
; _stencil_val(u.z,0,-1,0);
_stencil_val(ua.y,0,-1,1); _stencil_val(u.y,0,0,1);
_stencil_val(ua.y,0,-1,-1); _stencil_val(u.y,0,0,-1);


;;;

;;


;;
#line 167
_stencil_val_a(w.z,0,0,0);     

          
         
       
          
             
           


          
           
         
          
             
           

             

           


           

          
      }
  }end_foreach_level_or_leaf_stencil()

   BEGIN_FOREACH{
#line 163
foreach_level_or_leaf (l)

  {
    
      val(w.x,0,0,0) = (val(r.x,0,0,0)*sq(Delta) + dt/_const_rho*(2.*_const_mu.x*val(u.x,1,0,0) + 2.*_const_mu.x*val(u.x,-1,0,0)

        + _const_mu.y*(val(u.x,0,1,0) +
       (val(u.y,1,0,0) + val(ua.y,1,1,0))/4. -
       (val(u.y,-1,0,0) + val(ua.y,-1,1,0))/4.)
        - _const_mu.y*(- val(u.x,0,-1,0) +
           (val(ua.y,1,-1,0) + val(u.y,1,0,0))/4. -
           (val(ua.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


        + _const_mu.z*(val(u.x,0,0,1) +
         (val(u.z,1,0,0) + val(ua.z,1,0,1))/4. -
         (val(u.z,-1,0,0) + val(ua.z,-1,0,1))/4.)
        - _const_mu.z*(- val(u.x,0,0,-1) +
           (val(ua.z,1,0,-1) + val(u.z,1,0,0))/4. -
           (val(ua.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

        ))/
      (((coord){1.,1.,1.}).x*sq(Delta) + dt/_const_rho*(2.*_const_mu.x + 2.*_const_mu.x

          + _const_mu.y + _const_mu.y


          + _const_mu.z + _const_mu.z

          ));
      
#line 167
val(w.y,0,0,0) = (val(r.y,0,0,0)*sq(Delta) + dt/_const_rho*(2.*_const_mu.y*val(u.y,0,1,0) + 2.*_const_mu.y*val(u.y,0,-1,0)

        + _const_mu.z*(val(u.y,0,0,1) +
       (val(u.z,0,1,0) + val(ua.z,0,1,1))/4. -
       (val(u.z,0,-1,0) + val(ua.z,0,-1,1))/4.)
        - _const_mu.z*(- val(u.y,0,0,-1) +
           (val(ua.z,0,1,-1) + val(u.z,0,1,0))/4. -
           (val(ua.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


        + _const_mu.x*(val(u.y,1,0,0) +
         (val(u.x,0,1,0) + val(ua.x,1,1,0))/4. -
         (val(u.x,0,-1,0) + val(ua.x,1,-1,0))/4.)
        - _const_mu.x*(- val(u.y,-1,0,0) +
           (val(ua.x,-1,1,0) + val(u.x,0,1,0))/4. -
           (val(ua.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

        ))/
      (((coord){1.,1.,1.}).y*sq(Delta) + dt/_const_rho*(2.*_const_mu.y + 2.*_const_mu.y

          + _const_mu.z + _const_mu.z


          + _const_mu.x + _const_mu.x

          ));
      
#line 167
val(w.z,0,0,0) = (val(r.z,0,0,0)*sq(Delta) + dt/_const_rho*(2.*_const_mu.z*val(u.z,0,0,1) + 2.*_const_mu.z*val(u.z,0,0,-1)

        + _const_mu.x*(val(u.z,1,0,0) +
       (val(u.x,0,0,1) + val(ua.x,1,0,1))/4. -
       (val(u.x,0,0,-1) + val(ua.x,1,0,-1))/4.)
        - _const_mu.x*(- val(u.z,-1,0,0) +
           (val(ua.x,-1,0,1) + val(u.x,0,0,1))/4. -
           (val(ua.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


        + _const_mu.y*(val(u.z,0,1,0) +
         (val(u.y,0,0,1) + val(ua.y,0,1,1))/4. -
         (val(u.y,0,0,-1) + val(ua.y,0,1,-1))/4.)
        - _const_mu.y*(- val(u.z,0,-1,0) +
           (val(ua.y,0,-1,1) + val(u.y,0,0,1))/4. -
           (val(ua.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

        ))/
      (((coord){1.,1.,1.}).z*sq(Delta) + dt/_const_rho*(2.*_const_mu.z + 2.*_const_mu.z

          + _const_mu.x + _const_mu.x


          + _const_mu.y + _const_mu.y

          ));
  }end_foreach_level_or_leaf();}END_FOREACH }
#line 211 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
}
#line 220 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
static double residual_viscosity (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
          vector mu = p->mu;
          scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;
#line 266 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
  if(!is_constant(rho) && !is_constant(mu.x)){
#line 266 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
  foreach_stencil (1,{0},)
    { { 
_stencil_val(r.x,0,0,0);_stencil_val(u.x,0,0,0);
_stencil_val(rho,0,0,0);_stencil_val(mu.x,1,0,0);_stencil_val(u.x,1,0,0); _stencil_val(u.x,0,0,0);
_stencil_val(mu.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,-1,0,0); 

_stencil_val(mu.y,0,1,0);_stencil_val(u.x,0,1,0); _stencil_val(u.x,0,0,0);
_stencil_val(u.y,1,0,0); _stencil_val(u.y,1,1,0);
_stencil_val(u.y,-1,0,0); _stencil_val(u.y,-1,1,0); 
_stencil_val(mu.y,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,-1,0);
_stencil_val(u.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(u.y,-1,-1,0); _stencil_val(u.y,-1,0,0); 


_stencil_val(mu.z,0,0,1);_stencil_val(u.x,0,0,1); _stencil_val(u.x,0,0,0);
_stencil_val(u.z,1,0,0); _stencil_val(u.z,1,0,1);
_stencil_val(u.z,-1,0,0); _stencil_val(u.z,-1,0,1); 
_stencil_val(mu.z,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,0,-1);
_stencil_val(u.z,1,0,-1); _stencil_val(u.z,1,0,0);
_stencil_val(u.z,-1,0,-1); _stencil_val(u.z,-1,0,0);
      
#line 268
_stencil_val_a(res.x,0,0,0);
#line 288
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
} 
#line 267
{ 
_stencil_val(r.y,0,0,0);_stencil_val(u.y,0,0,0);
_stencil_val(rho,0,0,0);_stencil_val(mu.y,0,1,0);_stencil_val(u.y,0,1,0); _stencil_val(u.y,0,0,0);
_stencil_val(mu.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,-1,0); 

_stencil_val(mu.z,0,0,1);_stencil_val(u.y,0,0,1); _stencil_val(u.y,0,0,0);
_stencil_val(u.z,0,1,0); _stencil_val(u.z,0,1,1);
_stencil_val(u.z,0,-1,0); _stencil_val(u.z,0,-1,1); 
_stencil_val(mu.z,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0,-1);
_stencil_val(u.z,0,1,-1); _stencil_val(u.z,0,1,0);
_stencil_val(u.z,0,-1,-1); _stencil_val(u.z,0,-1,0); 


_stencil_val(mu.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(u.y,0,0,0);
_stencil_val(u.x,0,1,0); _stencil_val(u.x,1,1,0);
_stencil_val(u.x,0,-1,0); _stencil_val(u.x,1,-1,0); 
_stencil_val(mu.x,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,-1,0,0);
_stencil_val(u.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(u.x,-1,-1,0); _stencil_val(u.x,0,-1,0);
      
#line 268
_stencil_val_a(res.y,0,0,0);
#line 288
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
} 
#line 267
{ 
_stencil_val(r.z,0,0,0);_stencil_val(u.z,0,0,0);
_stencil_val(rho,0,0,0);_stencil_val(mu.z,0,0,1);_stencil_val(u.z,0,0,1); _stencil_val(u.z,0,0,0);
_stencil_val(mu.z,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,-1); 

_stencil_val(mu.x,1,0,0);_stencil_val(u.z,1,0,0); _stencil_val(u.z,0,0,0);
_stencil_val(u.x,0,0,1); _stencil_val(u.x,1,0,1);
_stencil_val(u.x,0,0,-1); _stencil_val(u.x,1,0,-1); 
_stencil_val(mu.x,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,-1,0,0);
_stencil_val(u.x,-1,0,1); _stencil_val(u.x,0,0,1);
_stencil_val(u.x,-1,0,-1); _stencil_val(u.x,0,0,-1); 


_stencil_val(mu.y,0,1,0);_stencil_val(u.z,0,1,0); _stencil_val(u.z,0,0,0);
_stencil_val(u.y,0,0,1); _stencil_val(u.y,0,1,1);
_stencil_val(u.y,0,0,-1); _stencil_val(u.y,0,1,-1); 
_stencil_val(mu.y,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,-1,0);
_stencil_val(u.y,0,-1,1); _stencil_val(u.y,0,0,1);
_stencil_val(u.y,0,-1,-1); _stencil_val(u.y,0,0,-1);
      
#line 268
_stencil_val_a(res.z,0,0,0);
#line 288
_stencil_val(res.z,0,0,0);
 {_stencil_val(res.z,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
}}end_foreach_stencil()
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(max:maxres)){
#line 266
foreach ()
    { {
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) +
        dt/val(rho,0,0,0)*(2.*val(mu.x,1,0,0)*(val(u.x,1,0,0) - val(u.x,0,0,0))
    - 2.*val(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))

    + val(mu.y,0,1,0)*(val(u.x,0,1,0) - val(u.x,0,0,0) +
          (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
          (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
    - val(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
       (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
       (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


    + val(mu.z,0,0,1)*(val(u.x,0,0,1) - val(u.x,0,0,0) +
     (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
     (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
    - val(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
       (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
       (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    } 
#line 267
{
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) +
        dt/val(rho,0,0,0)*(2.*val(mu.y,0,1,0)*(val(u.y,0,1,0) - val(u.y,0,0,0))
    - 2.*val(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))

    + val(mu.z,0,0,1)*(val(u.y,0,0,1) - val(u.y,0,0,0) +
          (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
          (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
    - val(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
       (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
       (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


    + val(mu.x,1,0,0)*(val(u.y,1,0,0) - val(u.y,0,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
    - val(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
       (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
       (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    } 
#line 267
{
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) +
        dt/val(rho,0,0,0)*(2.*val(mu.z,0,0,1)*(val(u.z,0,0,1) - val(u.z,0,0,0))
    - 2.*val(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))

    + val(mu.x,1,0,0)*(val(u.z,1,0,0) - val(u.z,0,0,0) +
          (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
          (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
    - val(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
       (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
       (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


    + val(mu.y,0,1,0)*(val(u.z,0,1,0) - val(u.z,0,0,0) +
     (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
     (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
    - val(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
       (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
       (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

    )/sq(Delta);
      if (fabs (val(res.z,0,0,0)) > maxres)
 maxres = fabs (val(res.z,0,0,0));
    }}end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 290
}else if(is_constant(rho) && !is_constant(mu.x)){double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);
#line 266 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
  foreach_stencil (1,{0},)
    { { 
_stencil_val(r.x,0,0,0);_stencil_val(u.x,0,0,0);
;_stencil_val(mu.x,1,0,0);_stencil_val(u.x,1,0,0); _stencil_val(u.x,0,0,0);
_stencil_val(mu.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,-1,0,0); 

_stencil_val(mu.y,0,1,0);_stencil_val(u.x,0,1,0); _stencil_val(u.x,0,0,0);
_stencil_val(u.y,1,0,0); _stencil_val(u.y,1,1,0);
_stencil_val(u.y,-1,0,0); _stencil_val(u.y,-1,1,0); 
_stencil_val(mu.y,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,-1,0);
_stencil_val(u.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(u.y,-1,-1,0); _stencil_val(u.y,-1,0,0); 


_stencil_val(mu.z,0,0,1);_stencil_val(u.x,0,0,1); _stencil_val(u.x,0,0,0);
_stencil_val(u.z,1,0,0); _stencil_val(u.z,1,0,1);
_stencil_val(u.z,-1,0,0); _stencil_val(u.z,-1,0,1); 
_stencil_val(mu.z,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,0,-1);
_stencil_val(u.z,1,0,-1); _stencil_val(u.z,1,0,0);
_stencil_val(u.z,-1,0,-1); _stencil_val(u.z,-1,0,0);
      
#line 268
_stencil_val_a(res.x,0,0,0);
#line 288
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
} 
#line 267
{ 
_stencil_val(r.y,0,0,0);_stencil_val(u.y,0,0,0);
;_stencil_val(mu.y,0,1,0);_stencil_val(u.y,0,1,0); _stencil_val(u.y,0,0,0);
_stencil_val(mu.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,-1,0); 

_stencil_val(mu.z,0,0,1);_stencil_val(u.y,0,0,1); _stencil_val(u.y,0,0,0);
_stencil_val(u.z,0,1,0); _stencil_val(u.z,0,1,1);
_stencil_val(u.z,0,-1,0); _stencil_val(u.z,0,-1,1); 
_stencil_val(mu.z,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0,-1);
_stencil_val(u.z,0,1,-1); _stencil_val(u.z,0,1,0);
_stencil_val(u.z,0,-1,-1); _stencil_val(u.z,0,-1,0); 


_stencil_val(mu.x,1,0,0);_stencil_val(u.y,1,0,0); _stencil_val(u.y,0,0,0);
_stencil_val(u.x,0,1,0); _stencil_val(u.x,1,1,0);
_stencil_val(u.x,0,-1,0); _stencil_val(u.x,1,-1,0); 
_stencil_val(mu.x,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,-1,0,0);
_stencil_val(u.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(u.x,-1,-1,0); _stencil_val(u.x,0,-1,0);
      
#line 268
_stencil_val_a(res.y,0,0,0);
#line 288
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
} 
#line 267
{ 
_stencil_val(r.z,0,0,0);_stencil_val(u.z,0,0,0);
;_stencil_val(mu.z,0,0,1);_stencil_val(u.z,0,0,1); _stencil_val(u.z,0,0,0);
_stencil_val(mu.z,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,-1); 

_stencil_val(mu.x,1,0,0);_stencil_val(u.z,1,0,0); _stencil_val(u.z,0,0,0);
_stencil_val(u.x,0,0,1); _stencil_val(u.x,1,0,1);
_stencil_val(u.x,0,0,-1); _stencil_val(u.x,1,0,-1); 
_stencil_val(mu.x,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,-1,0,0);
_stencil_val(u.x,-1,0,1); _stencil_val(u.x,0,0,1);
_stencil_val(u.x,-1,0,-1); _stencil_val(u.x,0,0,-1); 


_stencil_val(mu.y,0,1,0);_stencil_val(u.z,0,1,0); _stencil_val(u.z,0,0,0);
_stencil_val(u.y,0,0,1); _stencil_val(u.y,0,1,1);
_stencil_val(u.y,0,0,-1); _stencil_val(u.y,0,1,-1); 
_stencil_val(mu.y,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,-1,0);
_stencil_val(u.y,0,-1,1); _stencil_val(u.y,0,0,1);
_stencil_val(u.y,0,-1,-1); _stencil_val(u.y,0,0,-1);
      
#line 268
_stencil_val_a(res.z,0,0,0);
#line 288
_stencil_val(res.z,0,0,0);
 {_stencil_val(res.z,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
}}end_foreach_stencil()
#line 266 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(max:maxres)){
#line 266
foreach ()
    { {
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) +
        dt/_const_rho*(2.*val(mu.x,1,0,0)*(val(u.x,1,0,0) - val(u.x,0,0,0))
    - 2.*val(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))

    + val(mu.y,0,1,0)*(val(u.x,0,1,0) - val(u.x,0,0,0) +
          (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
          (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
    - val(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
       (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
       (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


    + val(mu.z,0,0,1)*(val(u.x,0,0,1) - val(u.x,0,0,0) +
     (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
     (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
    - val(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
       (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
       (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    } 
#line 267
{
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) +
        dt/_const_rho*(2.*val(mu.y,0,1,0)*(val(u.y,0,1,0) - val(u.y,0,0,0))
    - 2.*val(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))

    + val(mu.z,0,0,1)*(val(u.y,0,0,1) - val(u.y,0,0,0) +
          (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
          (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
    - val(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
       (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
       (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


    + val(mu.x,1,0,0)*(val(u.y,1,0,0) - val(u.y,0,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
    - val(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
       (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
       (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    } 
#line 267
{
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) +
        dt/_const_rho*(2.*val(mu.z,0,0,1)*(val(u.z,0,0,1) - val(u.z,0,0,0))
    - 2.*val(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))

    + val(mu.x,1,0,0)*(val(u.z,1,0,0) - val(u.z,0,0,0) +
          (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
          (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
    - val(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
       (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
       (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


    + val(mu.y,0,1,0)*(val(u.z,0,1,0) - val(u.z,0,0,0) +
     (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
     (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
    - val(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
       (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
       (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

    )/sq(Delta);
      if (fabs (val(res.z,0,0,0)) > maxres)
 maxres = fabs (val(res.z,0,0,0));
    }}end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 290
}else if(!is_constant(rho) && is_constant(mu.x)){_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX]};NOT_UNUSED(_const_mu);
#line 266 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
  foreach_stencil (1,{0},)
    { { 
_stencil_val(r.x,0,0,0);_stencil_val(u.x,0,0,0);
_stencil_val(rho,0,0,0);;_stencil_val(u.x,1,0,0); _stencil_val(u.x,0,0,0);
;_stencil_val(u.x,0,0,0); _stencil_val(u.x,-1,0,0);

;_stencil_val(u.x,0,1,0); _stencil_val(u.x,0,0,0);
_stencil_val(u.y,1,0,0); _stencil_val(u.y,1,1,0);
_stencil_val(u.y,-1,0,0); _stencil_val(u.y,-1,1,0);
;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,-1,0);
_stencil_val(u.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(u.y,-1,-1,0); _stencil_val(u.y,-1,0,0);


;_stencil_val(u.x,0,0,1); _stencil_val(u.x,0,0,0);
_stencil_val(u.z,1,0,0); _stencil_val(u.z,1,0,1);
_stencil_val(u.z,-1,0,0); _stencil_val(u.z,-1,0,1);
;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,0,-1);
_stencil_val(u.z,1,0,-1); _stencil_val(u.z,1,0,0);
_stencil_val(u.z,-1,0,-1); _stencil_val(u.z,-1,0,0);
      
#line 268
_stencil_val_a(res.x,0,0,0);
#line 288
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
} 
#line 267
{ 
_stencil_val(r.y,0,0,0);_stencil_val(u.y,0,0,0);
_stencil_val(rho,0,0,0);;_stencil_val(u.y,0,1,0); _stencil_val(u.y,0,0,0);
;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,-1,0);

;_stencil_val(u.y,0,0,1); _stencil_val(u.y,0,0,0);
_stencil_val(u.z,0,1,0); _stencil_val(u.z,0,1,1);
_stencil_val(u.z,0,-1,0); _stencil_val(u.z,0,-1,1);
;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0,-1);
_stencil_val(u.z,0,1,-1); _stencil_val(u.z,0,1,0);
_stencil_val(u.z,0,-1,-1); _stencil_val(u.z,0,-1,0);


;_stencil_val(u.y,1,0,0); _stencil_val(u.y,0,0,0);
_stencil_val(u.x,0,1,0); _stencil_val(u.x,1,1,0);
_stencil_val(u.x,0,-1,0); _stencil_val(u.x,1,-1,0);
;_stencil_val(u.y,0,0,0); _stencil_val(u.y,-1,0,0);
_stencil_val(u.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(u.x,-1,-1,0); _stencil_val(u.x,0,-1,0);
      
#line 268
_stencil_val_a(res.y,0,0,0);
#line 288
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
} 
#line 267
{ 
_stencil_val(r.z,0,0,0);_stencil_val(u.z,0,0,0);
_stencil_val(rho,0,0,0);;_stencil_val(u.z,0,0,1); _stencil_val(u.z,0,0,0);
;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,-1);

;_stencil_val(u.z,1,0,0); _stencil_val(u.z,0,0,0);
_stencil_val(u.x,0,0,1); _stencil_val(u.x,1,0,1);
_stencil_val(u.x,0,0,-1); _stencil_val(u.x,1,0,-1);
;_stencil_val(u.z,0,0,0); _stencil_val(u.z,-1,0,0);
_stencil_val(u.x,-1,0,1); _stencil_val(u.x,0,0,1);
_stencil_val(u.x,-1,0,-1); _stencil_val(u.x,0,0,-1);


;_stencil_val(u.z,0,1,0); _stencil_val(u.z,0,0,0);
_stencil_val(u.y,0,0,1); _stencil_val(u.y,0,1,1);
_stencil_val(u.y,0,0,-1); _stencil_val(u.y,0,1,-1);
;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,-1,0);
_stencil_val(u.y,0,-1,1); _stencil_val(u.y,0,0,1);
_stencil_val(u.y,0,-1,-1); _stencil_val(u.y,0,0,-1);
      
#line 268
_stencil_val_a(res.z,0,0,0);
#line 288
_stencil_val(res.z,0,0,0);
 {_stencil_val(res.z,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
}}end_foreach_stencil()
#line 266 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(max:maxres)){
#line 266
foreach ()
    { {
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) +
        dt/val(rho,0,0,0)*(2.*_const_mu.x*(val(u.x,1,0,0) - val(u.x,0,0,0))
    - 2.*_const_mu.x*(val(u.x,0,0,0) - val(u.x,-1,0,0))

    + _const_mu.y*(val(u.x,0,1,0) - val(u.x,0,0,0) +
          (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
          (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
    - _const_mu.y*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
       (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
       (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


    + _const_mu.z*(val(u.x,0,0,1) - val(u.x,0,0,0) +
     (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
     (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
    - _const_mu.z*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
       (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
       (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    } 
#line 267
{
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) +
        dt/val(rho,0,0,0)*(2.*_const_mu.y*(val(u.y,0,1,0) - val(u.y,0,0,0))
    - 2.*_const_mu.y*(val(u.y,0,0,0) - val(u.y,0,-1,0))

    + _const_mu.z*(val(u.y,0,0,1) - val(u.y,0,0,0) +
          (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
          (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
    - _const_mu.z*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
       (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
       (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


    + _const_mu.x*(val(u.y,1,0,0) - val(u.y,0,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
    - _const_mu.x*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
       (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
       (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    } 
#line 267
{
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) +
        dt/val(rho,0,0,0)*(2.*_const_mu.z*(val(u.z,0,0,1) - val(u.z,0,0,0))
    - 2.*_const_mu.z*(val(u.z,0,0,0) - val(u.z,0,0,-1))

    + _const_mu.x*(val(u.z,1,0,0) - val(u.z,0,0,0) +
          (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
          (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
    - _const_mu.x*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
       (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
       (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


    + _const_mu.y*(val(u.z,0,1,0) - val(u.z,0,0,0) +
     (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
     (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
    - _const_mu.y*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
       (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
       (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

    )/sq(Delta);
      if (fabs (val(res.z,0,0,0)) > maxres)
 maxres = fabs (val(res.z,0,0,0));
    }}end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 290
}else {double _const_rho=_constant[rho.i-_NVARMAX];NOT_UNUSED(_const_rho);_coord _const_mu={_constant[mu.x.i-_NVARMAX],_constant[mu.y.i-_NVARMAX],_constant[mu.z.i-_NVARMAX]};NOT_UNUSED(_const_mu);
#line 266 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
  foreach_stencil (1,{0},)
    { { 
_stencil_val(r.x,0,0,0);_stencil_val(u.x,0,0,0);
;;_stencil_val(u.x,1,0,0); _stencil_val(u.x,0,0,0);
;_stencil_val(u.x,0,0,0); _stencil_val(u.x,-1,0,0);

;_stencil_val(u.x,0,1,0); _stencil_val(u.x,0,0,0);
_stencil_val(u.y,1,0,0); _stencil_val(u.y,1,1,0);
_stencil_val(u.y,-1,0,0); _stencil_val(u.y,-1,1,0);
;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,-1,0);
_stencil_val(u.y,1,-1,0); _stencil_val(u.y,1,0,0);
_stencil_val(u.y,-1,-1,0); _stencil_val(u.y,-1,0,0);


;_stencil_val(u.x,0,0,1); _stencil_val(u.x,0,0,0);
_stencil_val(u.z,1,0,0); _stencil_val(u.z,1,0,1);
_stencil_val(u.z,-1,0,0); _stencil_val(u.z,-1,0,1);
;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0,0,-1);
_stencil_val(u.z,1,0,-1); _stencil_val(u.z,1,0,0);
_stencil_val(u.z,-1,0,-1); _stencil_val(u.z,-1,0,0);
      
#line 268
_stencil_val_a(res.x,0,0,0);
#line 288
_stencil_val(res.x,0,0,0);
 {_stencil_val(res.x,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
} 
#line 267
{ 
_stencil_val(r.y,0,0,0);_stencil_val(u.y,0,0,0);
;;_stencil_val(u.y,0,1,0); _stencil_val(u.y,0,0,0);
;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,-1,0);

;_stencil_val(u.y,0,0,1); _stencil_val(u.y,0,0,0);
_stencil_val(u.z,0,1,0); _stencil_val(u.z,0,1,1);
_stencil_val(u.z,0,-1,0); _stencil_val(u.z,0,-1,1);
;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0,-1);
_stencil_val(u.z,0,1,-1); _stencil_val(u.z,0,1,0);
_stencil_val(u.z,0,-1,-1); _stencil_val(u.z,0,-1,0);


;_stencil_val(u.y,1,0,0); _stencil_val(u.y,0,0,0);
_stencil_val(u.x,0,1,0); _stencil_val(u.x,1,1,0);
_stencil_val(u.x,0,-1,0); _stencil_val(u.x,1,-1,0);
;_stencil_val(u.y,0,0,0); _stencil_val(u.y,-1,0,0);
_stencil_val(u.x,-1,1,0); _stencil_val(u.x,0,1,0);
_stencil_val(u.x,-1,-1,0); _stencil_val(u.x,0,-1,0);
      
#line 268
_stencil_val_a(res.y,0,0,0);
#line 288
_stencil_val(res.y,0,0,0);
 {_stencil_val(res.y,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
} 
#line 267
{ 
_stencil_val(r.z,0,0,0);_stencil_val(u.z,0,0,0);
;;_stencil_val(u.z,0,0,1); _stencil_val(u.z,0,0,0);
;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,-1);

;_stencil_val(u.z,1,0,0); _stencil_val(u.z,0,0,0);
_stencil_val(u.x,0,0,1); _stencil_val(u.x,1,0,1);
_stencil_val(u.x,0,0,-1); _stencil_val(u.x,1,0,-1);
;_stencil_val(u.z,0,0,0); _stencil_val(u.z,-1,0,0);
_stencil_val(u.x,-1,0,1); _stencil_val(u.x,0,0,1);
_stencil_val(u.x,-1,0,-1); _stencil_val(u.x,0,0,-1);


;_stencil_val(u.z,0,1,0); _stencil_val(u.z,0,0,0);
_stencil_val(u.y,0,0,1); _stencil_val(u.y,0,1,1);
_stencil_val(u.y,0,0,-1); _stencil_val(u.y,0,1,-1);
;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,-1,0);
_stencil_val(u.y,0,-1,1); _stencil_val(u.y,0,0,1);
_stencil_val(u.y,0,-1,-1); _stencil_val(u.y,0,0,-1);
      
#line 268
_stencil_val_a(res.z,0,0,0);
#line 288
_stencil_val(res.z,0,0,0);
 {_stencil_val(res.z,0,0,0);   }    
         
      

       
            
          
       
         
       


       
       
     
       
         
       

    
          
    
#line 290
}}end_foreach_stencil()
#line 266 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(max:maxres)){
#line 266
foreach ()
    { {
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) +
        dt/_const_rho*(2.*_const_mu.x*(val(u.x,1,0,0) - val(u.x,0,0,0))
    - 2.*_const_mu.x*(val(u.x,0,0,0) - val(u.x,-1,0,0))

    + _const_mu.y*(val(u.x,0,1,0) - val(u.x,0,0,0) +
          (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
          (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
    - _const_mu.y*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
       (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
       (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


    + _const_mu.z*(val(u.x,0,0,1) - val(u.x,0,0,0) +
     (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
     (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
    - _const_mu.z*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
       (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
       (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.x,0,0,0)) > maxres)
 maxres = fabs (val(res.x,0,0,0));
    } 
#line 267
{
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) +
        dt/_const_rho*(2.*_const_mu.y*(val(u.y,0,1,0) - val(u.y,0,0,0))
    - 2.*_const_mu.y*(val(u.y,0,0,0) - val(u.y,0,-1,0))

    + _const_mu.z*(val(u.y,0,0,1) - val(u.y,0,0,0) +
          (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
          (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
    - _const_mu.z*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
       (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
       (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


    + _const_mu.x*(val(u.y,1,0,0) - val(u.y,0,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
    - _const_mu.x*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
       (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
       (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.y,0,0,0)) > maxres)
 maxres = fabs (val(res.y,0,0,0));
    } 
#line 267
{
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) +
        dt/_const_rho*(2.*_const_mu.z*(val(u.z,0,0,1) - val(u.z,0,0,0))
    - 2.*_const_mu.z*(val(u.z,0,0,0) - val(u.z,0,0,-1))

    + _const_mu.x*(val(u.z,1,0,0) - val(u.z,0,0,0) +
          (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
          (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
    - _const_mu.x*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
       (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
       (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


    + _const_mu.y*(val(u.z,0,1,0) - val(u.z,0,0,0) +
     (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
     (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
    - _const_mu.y*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
       (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
       (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

    )/sq(Delta);
      if (fabs (val(res.z,0,0,0)) > maxres)
 maxres = fabs (val(res.z,0,0,0));
    }}end_foreach();mpi_all_reduce_array(&maxres,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 290
}

  return maxres;
}
#line 304 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
     
mgstats viscosity (vector u, vector mu, scalar rho, double dt,
     int nrelax, scalar * res)
{tracing("viscosity","/home/esteban/ProgramFile/basilisk/src/viscosity.h",305);





  vector  r=new_vector("r");
  foreach_stencil(1,{0},)
    {
      { _stencil_val(u.x,0,0,0);_stencil_val_a(r.x,0,0,0); }
      
#line 316
{ _stencil_val(u.y,0,0,0);_stencil_val_a(r.y,0,0,0); }
      
#line 316
{ _stencil_val(u.z,0,0,0);_stencil_val_a(r.z,0,0,0); }}end_foreach_stencil()
   BEGIN_FOREACH{
#line 314
foreach()
    {
      val(r.x,0,0,0) = val(u.x,0,0,0);
      
#line 316
val(r.y,0,0,0) = val(u.y,0,0,0);
      
#line 316
val(r.z,0,0,0) = val(u.z,0,0,0);}end_foreach();}END_FOREACH 




  restriction (((scalar[]){mu.x,mu.y,mu.z,rho,{-1}}));
  struct Viscosity p = { mu, rho, dt };
  { mgstats _ret= mg_solve ((scalar *)((vector[]){u,{{-1},{-1},{-1}}}), (scalar *)((vector[]){r,{{-1},{-1},{-1}}})
,
     
#line 324
residual_viscosity, relax_viscosity, &p, nrelax, res
#line 135 "/home/esteban/ProgramFile/basilisk/src/poisson.h"
, 
0, 
TOLERANCE
#line 324 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
);delete((scalar*)((vector[]){r,{{-1},{-1},{-1}}}));{end_tracing("viscosity","/home/esteban/ProgramFile/basilisk/src/viscosity.h",324);return _ret;}}delete((scalar*)((vector[]){r,{{-1},{-1},{-1}}}));
end_tracing("viscosity","/home/esteban/ProgramFile/basilisk/src/viscosity.h",325);}
#line 334 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
     
mgstats viscosity_explicit (vector u, vector mu, scalar rho, double dt)
{tracing("viscosity_explicit","/home/esteban/ProgramFile/basilisk/src/viscosity.h",335);
  vector  r=new_vector("r");
  mgstats mg = {0};
  struct Viscosity p = { mu, rho, dt };
  mg.resb = residual_viscosity ((scalar *)((vector[]){u,{{-1},{-1},{-1}}}), (scalar *)((vector[]){u,{{-1},{-1},{-1}}}), (scalar *)((vector[]){r,{{-1},{-1},{-1}}}), &p);
  foreach_stencil(1,{0},)
    {
      { _stencil_val(r.x,0,0,0);_stencil_val_r(u.x,0,0,0); }
      
#line 343
{ _stencil_val(r.y,0,0,0);_stencil_val_r(u.y,0,0,0); }
      
#line 343
{ _stencil_val(r.z,0,0,0);_stencil_val_r(u.z,0,0,0); }}end_foreach_stencil()
   BEGIN_FOREACH{
#line 341
foreach()
    {
      val(u.x,0,0,0) += val(r.x,0,0,0);
      
#line 343
val(u.y,0,0,0) += val(r.y,0,0,0);
      
#line 343
val(u.z,0,0,0) += val(r.z,0,0,0);}end_foreach();}END_FOREACH 
  {delete((scalar*)((vector[]){r,{{-1},{-1},{-1}}}));{end_tracing("viscosity_explicit","/home/esteban/ProgramFile/basilisk/src/viscosity.h",344);return mg;}}delete((scalar*)((vector[]){r,{{-1},{-1},{-1}}}));
end_tracing("viscosity_explicit","/home/esteban/ProgramFile/basilisk/src/viscosity.h",345);}
#line 34 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
#line 44 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
scalar  p={0};
vector  u={{1},{2},{3}},  g={{4},{5},{6}};
scalar  pf={7};
vector  uf={{8},{9},{10}};
#line 70 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
        vector mu = {{_NVARMAX+4},{_NVARMAX+5},{_NVARMAX+6}}, a = {{_NVARMAX+4},{_NVARMAX+5},{_NVARMAX+6}}, alpha = {{_NVARMAX+7},{_NVARMAX+8},{_NVARMAX+9}};
        scalar rho = {_NVARMAX+10};
mgstats mgp = {0}, mgpf = {0}, mgu = {0};
bool stokes = false;
#line 91
static double _boundary0(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.x*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.x*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.x,1,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}static double _boundary0_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*val(fm.x,1,0,0)/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*_const_fm.x/val(alpha.x,1,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*val(fm.x,1,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.x,1,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}
static double _boundary1(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.x,0,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}static double _boundary1_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*val(fm.x,0,0,0)/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*_const_fm.x/val(alpha.x,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*val(fm.x,0,0,0)/_const_alpha.x), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.x,0,0,0)*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.x*_const_fm.x/_const_alpha.x), point, neighbor, _s, data));}}}}}








static double _boundary2(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.y*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.y*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.y,0,1,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}static double _boundary2_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*val(fm.y,0,1,0)/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*_const_fm.y/val(alpha.y,0,1,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*val(fm.y,0,1,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.y,0,1,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}
static double _boundary3(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.y,0,0,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}static double _boundary3_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*val(fm.y,0,0,0)/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*_const_fm.y/val(alpha.y,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*val(fm.y,0,0,0)/_const_alpha.y), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.y,0,0,0)*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.y*_const_fm.y/_const_alpha.y), point, neighbor, _s, data));}}}}}


static double _boundary4(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.z,0,0,1)*val(fm.z,0,0,1)/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.z*val(fm.z,0,0,1)/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.z,0,0,1)*_const_fm.z/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.z*_const_fm.z/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.z,0,0,1)*val(fm.z,0,0,1)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.z*val(fm.z,0,0,1)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((val(a.z,0,0,1)*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann((_const_a.z*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}}static double _boundary4_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.z,0,0,1)*val(fm.z,0,0,1)/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.z*val(fm.z,0,0,1)/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.z,0,0,1)*_const_fm.z/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.z*_const_fm.z/val(alpha.z,0,0,1)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.z,0,0,1)*val(fm.z,0,0,1)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.z*val(fm.z,0,0,1)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((val(a.z,0,0,1)*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous((_const_a.z*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}}
static double _boundary5(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.z,0,0,0)*val(fm.z,0,0,0)/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.z*val(fm.z,0,0,0)/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.z,0,0,0)*_const_fm.z/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.z*_const_fm.z/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.z,0,0,0)*val(fm.z,0,0,0)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.z*val(fm.z,0,0,0)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (val(a.z,0,0,0)*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(- (_const_a.z*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}}static double _boundary5_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;if(!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.z,0,0,0)*val(fm.z,0,0,0)/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.z*val(fm.z,0,0,0)/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.z,0,0,0)*_const_fm.z/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.z*_const_fm.z/val(alpha.z,0,0,0)), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.z,0,0,0)*val(fm.z,0,0,0)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.z*val(fm.z,0,0,0)/_const_alpha.z), point, neighbor, _s, data));}}}}else if(!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (val(a.z,0,0,0)*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(- (_const_a.z*_const_fm.z/_const_alpha.z), point, neighbor, _s, data));}}}}}
#line 126
static int defaults_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}
#line 126 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
      static int defaults_0(const int i,const double t,Event *_ev){tracing("defaults_0","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",126);
{




  mgp = (mgstats){0};
  mgpf = (mgstats){0};
  mgu = (mgstats){0};

  CFL = 0.8;




  _attribute[p.i].nodump = _attribute[pf.i].nodump = true;





  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    vector alphav = alpha;
    if(!is_constant(fm.x)){
    
#line 153
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
      { _stencil_val(fm.x,0,0,0);_stencil_val_a(alphav.x,0,0,0); }}end__stencil_is_face_x()
#line 153
_stencil_is_face_y(){
      { _stencil_val(fm.y,0,0,0);_stencil_val_a(alphav.y,0,0,0); }}end__stencil_is_face_y()
#line 153
_stencil_is_face_z(){
      { _stencil_val(fm.z,0,0,0);_stencil_val_a(alphav.z,0,0,0); }}end__stencil_is_face_z()}end_foreach_face_stencil() BEGIN_FOREACH{
#line 153
foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = val(fm.x,0,0,0);}end_is_face_x()
#line 153
is_face_y(){
      val(alphav.y,0,0,0) = val(fm.y,0,0,0);}end_is_face_y()
#line 153
is_face_z(){
      val(alphav.z,0,0,0) = val(fm.z,0,0,0);}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
    
#line 153
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
      {;_stencil_val_a(alphav.x,0,0,0); }}end__stencil_is_face_x()
#line 153
_stencil_is_face_y(){
      {;_stencil_val_a(alphav.y,0,0,0); }}end__stencil_is_face_y()
#line 153
_stencil_is_face_z(){
      {;_stencil_val_a(alphav.z,0,0,0); }}end__stencil_is_face_z()}end_foreach_face_stencil()
     BEGIN_FOREACH{
#line 153
foreach_face_generic(){is_face_x(){
      val(alphav.x,0,0,0) = _const_fm.x;}end_is_face_x()
#line 153
is_face_y(){
      val(alphav.y,0,0,0) = _const_fm.y;}end_is_face_y()
#line 153
is_face_z(){
      val(alphav.z,0,0,0) = _const_fm.z;}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }
  }
#line 185 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
  foreach_stencil(1,{0},)
    {
      {_stencil_val(u.x,0,0,0);   }
      
#line 187
{_stencil_val(u.y,0,0,0);   }
      
#line 187
{_stencil_val(u.z,0,0,0);   }}end_foreach_stencil()
#line 185 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
   BEGIN_FOREACH{foreach()
    {
      dimensional (val(u.x,0,0,0) == Delta/t);
      
#line 187
dimensional (val(u.y,0,0,0) == Delta/t);
      
#line 187
dimensional (val(u.z,0,0,0) == Delta/t);}end_foreach();}END_FOREACH 
}{end_tracing("defaults_0","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",188);return 0;}end_tracing("defaults_0","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",188);}





static int default_display_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}






#line 194
      static int default_display(const int i,const double t,Event *_ev){tracing("default_display","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",194);
{display ("squares (color = 'u.x', spread = -1);"
#line 494 "/home/esteban/ProgramFile/basilisk/src/common.h"
, false
#line 195 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
);
  
#line 195
}{end_tracing("default_display","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",195);return 0;}end_tracing("default_display","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",195);}





double dtmax;

static int init_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}


#line 203
      static int init(const int i,const double t,Event *_ev){tracing("init","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",203);
{
  trash (((vector[]){uf,{{-1},{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 206
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val_a(uf.x,0,0,0);  }}end__stencil_is_face_x()
#line 206
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val_a(uf.y,0,0,0);  }}end__stencil_is_face_y()
#line 206
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);_stencil_val_a(uf.z,0,0,0);  }}end__stencil_is_face_z()}end_foreach_face_stencil() BEGIN_FOREACH{
#line 206
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()
#line 206
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()
#line 206
is_face_z(){
    val(uf.z,0,0,0) = val(fm.z,0,0,0)*((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.);}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 206
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val_a(uf.x,0,0,0);  }}end__stencil_is_face_x()
#line 206
_stencil_is_face_y(){
    {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val_a(uf.y,0,0,0);  }}end__stencil_is_face_y()
#line 206
_stencil_is_face_z(){
    {;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);_stencil_val_a(uf.z,0,0,0);  }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 206
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.);}end_is_face_x()
#line 206
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.);}end_is_face_y()
#line 206
is_face_z(){
    val(uf.z,0,0,0) = _const_fm.z*((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.);}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }




  event ("properties");





  dtmax = DT;
  event ("stability");
}{end_tracing("init","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",220);return 0;}end_tracing("init","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",220);}








static int set_dtmax_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 229 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
      static int set_dtmax(const int i,const double t,Event *_ev){tracing("set_dtmax","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",229);{dtmax = DT; }{end_tracing("set_dtmax","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",229);return 0;}end_tracing("set_dtmax","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",229);}

static int stability_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 231
      static int stability(const int i,const double t,Event *_ev){tracing("stability","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",231); {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
}{end_tracing("stability","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",233);return 0;}end_tracing("stability","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",233);}







static int vof_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}








#line 241
static int vof(const int i,const double t,Event *_ev){;return 0;}
static int tracer_advection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}

#line 242
static int tracer_advection(const int i,const double t,Event *_ev){;return 0;}
static int tracer_diffusion_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}

#line 243
static int tracer_diffusion(const int i,const double t,Event *_ev){;return 0;}






static int properties_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}







#line 250
static int properties(const int i,const double t,Event *_ev){;return 0;}
#line 262 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
void prediction()
{
  vector du;
   {
    scalar s = new_scalar("s");
    du.x = s;
  } 
#line 265
{
    scalar s = new_scalar("s");
    du.y = s;
  } 
#line 265
{
    scalar s = new_scalar("s");
    du.z = s;
  }

  if (_attribute[u.x.i].gradient)
    {
    
#line 271
foreach_stencil(1,{0},)
      { {





_stencil_val(u.x,-1,0,0); _stencil_val(u.x,0,0,0); _stencil_val(u.x,1,0,0);





   
#line 278
_stencil_val_a(du.x,0,0,0);   
      } 
#line 272
{





_stencil_val(u.y,0,-1,0); _stencil_val(u.y,0,0,0); _stencil_val(u.y,0,1,0);





   
#line 278
_stencil_val_a(du.y,0,0,0);   
      } 
#line 272
{





_stencil_val(u.z,0,0,-1); _stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,1);





   
#line 278
_stencil_val_a(du.z,0,0,0);   
      }}end_foreach_stencil() BEGIN_FOREACH{
#line 271
foreach()
      { {





   val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
      } 
#line 272
{





   val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;
      } 
#line 272
{





   val(du.z,0,0,0) = _attribute[u.z.i].gradient (val(u.z,0,0,-1), val(u.z,0,0,0), val(u.z,0,0,1))/Delta;
      }}end_foreach();}END_FOREACH }
  else
    {
    
#line 281
foreach_stencil(1,{0},)
      { {





_stencil_val(u.x,1,0,0); _stencil_val(u.x,-1,0,0);





   
#line 288
_stencil_val_a(du.x,0,0,0);   
    } 
#line 282
{





_stencil_val(u.y,0,1,0); _stencil_val(u.y,0,-1,0);





   
#line 288
_stencil_val_a(du.y,0,0,0);   
    } 
#line 282
{





_stencil_val(u.z,0,0,1); _stencil_val(u.z,0,0,-1);





   
#line 288
_stencil_val_a(du.z,0,0,0);   
    }}end_foreach_stencil() BEGIN_FOREACH{
#line 281
foreach()
      { {





   val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
    } 
#line 282
{





   val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);
    } 
#line 282
{





   val(du.z,0,0,0) = (val(u.z,0,0,1) - val(u.z,0,0,-1))/(2.*Delta);
    }}end_foreach();}END_FOREACH }

  trash (((vector[]){uf,{{-1},{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 292
foreach_face_stencil(1,{0},){_stencil_is_face_x(){ {       
     _stencil_val(u.x,-1,0,0);_stencil_val(u.x,0,0,0);     
    
_stencil_val(u.x, o_stencil,0,0);_stencil_val(g.x,0,0,0); _stencil_val(g.x,-1,0,0);_stencil_val(du.x,o_stencil,0,0);
    
#line 295
_stencil_val_a(uf.x,0,0,0);

_stencil_val(fm.y,o_stencil,0,0); _stencil_val(fm.y,o_stencil,1,0); {        
       _stencil_val(u.x,o_stencil,-1,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,1,0);_stencil_val(u.y, o_stencil,0,0);
_stencil_val(u.y,o_stencil,0,0);
      
#line 299
_stencil_val_r(uf.x,0,0,0);  
    }


_stencil_val(fm.z,o_stencil,0,0); _stencil_val(fm.z,o_stencil,0,1); {        
       _stencil_val(u.x,o_stencil,0,-1);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,0,1);_stencil_val(u.z, o_stencil,0,0);
_stencil_val(u.z,o_stencil,0,0);
      
#line 305
_stencil_val_r(uf.x,0,0,0);  
    } 

_stencil_val(fm.x,0,0,0);        

      


      

    
#line 308
_stencil_val_r(uf.x,0,0,0); 
  }}end__stencil_is_face_x()
#line 292
_stencil_is_face_y(){ {       
     _stencil_val(u.y,0,-1,0);_stencil_val(u.y,0,0,0);     
    
_stencil_val(u.y,0, o_stencil,0);_stencil_val(g.y,0,0,0); _stencil_val(g.y,0,-1,0);_stencil_val(du.y,0,o_stencil,0);
    
#line 295
_stencil_val_a(uf.y,0,0,0);

_stencil_val(fm.z,0,o_stencil,0); _stencil_val(fm.z,0,o_stencil,1); {        
       _stencil_val(u.y,0,o_stencil,-1);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,0,o_stencil,1);_stencil_val(u.z,0, o_stencil,0);
_stencil_val(u.z,0,o_stencil,0);
      
#line 299
_stencil_val_r(uf.y,0,0,0);  
    }


_stencil_val(fm.x,0,o_stencil,0); _stencil_val(fm.x,1,o_stencil,0); {        
       _stencil_val(u.y,-1,o_stencil,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,1,o_stencil,0);_stencil_val(u.x,0, o_stencil,0);
_stencil_val(u.x,0,o_stencil,0);
      
#line 305
_stencil_val_r(uf.y,0,0,0);  
    } 

_stencil_val(fm.y,0,0,0);        

      


      

    
#line 308
_stencil_val_r(uf.y,0,0,0); 
  }}end__stencil_is_face_y()
#line 292
_stencil_is_face_z(){ {       
     _stencil_val(u.z,0,0,-1);_stencil_val(u.z,0,0,0);     
    
_stencil_val(u.z,0,0, o_stencil);_stencil_val(g.z,0,0,0); _stencil_val(g.z,0,0,-1);_stencil_val(du.z,0,0,o_stencil);
    
#line 295
_stencil_val_a(uf.z,0,0,0);

_stencil_val(fm.x,0,0,o_stencil); _stencil_val(fm.x,1,0,o_stencil); {        
       _stencil_val(u.z,-1,0,o_stencil);_stencil_val(u.z,0,0, o_stencil);_stencil_val(u.z,0,0, o_stencil); _stencil_val(u.z,1,0,o_stencil);_stencil_val(u.x,0,0, o_stencil);
_stencil_val(u.x,0,0,o_stencil);
      
#line 299
_stencil_val_r(uf.z,0,0,0);  
    }


_stencil_val(fm.y,0,0,o_stencil); _stencil_val(fm.y,0,1,o_stencil); {        
       _stencil_val(u.z,0,-1,o_stencil);_stencil_val(u.z,0,0, o_stencil);_stencil_val(u.z,0,0, o_stencil); _stencil_val(u.z,0,1,o_stencil);_stencil_val(u.y,0,0, o_stencil);
_stencil_val(u.y,0,0,o_stencil);
      
#line 305
_stencil_val_r(uf.z,0,0,0);  
    } 

_stencil_val(fm.z,0,0,0);        

      


      

    
#line 308
_stencil_val_r(uf.z,0,0,0); 
  }}end__stencil_is_face_z()}end_foreach_face_stencil() BEGIN_FOREACH{
#line 292
foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (val(fm.y,i,0,0) && val(fm.y,i,1,0)) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }


    if (val(fm.z,i,0,0) && val(fm.z,i,0,1)) {
      double fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);
    }

    val(uf.x,0,0,0) *= val(fm.x,0,0,0);
  }}end_is_face_x()
#line 292
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (val(fm.z,0,i,0) && val(fm.z,0,i,1)) {
      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);
    }


    if (val(fm.x,0,i,0) && val(fm.x,1,i,0)) {
      double fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);
    }

    val(uf.y,0,0,0) *= val(fm.y,0,0,0);
  }}end_is_face_y()
#line 292
is_face_z(){ {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;

    if (val(fm.x,0,0,i) && val(fm.x,1,0,i)) {
      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);
    }


    if (val(fm.y,0,0,i) && val(fm.y,0,1,i)) {
      double fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);
    }

    val(uf.z,0,0,0) *= val(fm.z,0,0,0);
  }}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 292
foreach_face_stencil(1,{0},){_stencil_is_face_x(){ {       
     _stencil_val(u.x,-1,0,0);_stencil_val(u.x,0,0,0);     
    
_stencil_val(u.x, o_stencil,0,0);_stencil_val(g.x,0,0,0); _stencil_val(g.x,-1,0,0);_stencil_val(du.x,o_stencil,0,0);
    
#line 295
_stencil_val_a(uf.x,0,0,0);

;; {        
       _stencil_val(u.x,o_stencil,-1,0);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,1,0);_stencil_val(u.y, o_stencil,0,0);
_stencil_val(u.y,o_stencil,0,0);
      
#line 299
_stencil_val_r(uf.x,0,0,0);  
    }


;; {        
       _stencil_val(u.x,o_stencil,0,-1);_stencil_val(u.x, o_stencil,0,0);_stencil_val(u.x, o_stencil,0,0); _stencil_val(u.x,o_stencil,0,1);_stencil_val(u.z, o_stencil,0,0);
_stencil_val(u.z,o_stencil,0,0);
      
#line 305
_stencil_val_r(uf.x,0,0,0);  
    }

;        

      


      

    
#line 308
_stencil_val_r(uf.x,0,0,0); 
  }}end__stencil_is_face_x()
#line 292
_stencil_is_face_y(){ {       
     _stencil_val(u.y,0,-1,0);_stencil_val(u.y,0,0,0);     
    
_stencil_val(u.y,0, o_stencil,0);_stencil_val(g.y,0,0,0); _stencil_val(g.y,0,-1,0);_stencil_val(du.y,0,o_stencil,0);
    
#line 295
_stencil_val_a(uf.y,0,0,0);

;; {        
       _stencil_val(u.y,0,o_stencil,-1);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,0,o_stencil,1);_stencil_val(u.z,0, o_stencil,0);
_stencil_val(u.z,0,o_stencil,0);
      
#line 299
_stencil_val_r(uf.y,0,0,0);  
    }


;; {        
       _stencil_val(u.y,-1,o_stencil,0);_stencil_val(u.y,0, o_stencil,0);_stencil_val(u.y,0, o_stencil,0); _stencil_val(u.y,1,o_stencil,0);_stencil_val(u.x,0, o_stencil,0);
_stencil_val(u.x,0,o_stencil,0);
      
#line 305
_stencil_val_r(uf.y,0,0,0);  
    }

;        

      


      

    
#line 308
_stencil_val_r(uf.y,0,0,0); 
  }}end__stencil_is_face_y()
#line 292
_stencil_is_face_z(){ {       
     _stencil_val(u.z,0,0,-1);_stencil_val(u.z,0,0,0);     
    
_stencil_val(u.z,0,0, o_stencil);_stencil_val(g.z,0,0,0); _stencil_val(g.z,0,0,-1);_stencil_val(du.z,0,0,o_stencil);
    
#line 295
_stencil_val_a(uf.z,0,0,0);

;; {        
       _stencil_val(u.z,-1,0,o_stencil);_stencil_val(u.z,0,0, o_stencil);_stencil_val(u.z,0,0, o_stencil); _stencil_val(u.z,1,0,o_stencil);_stencil_val(u.x,0,0, o_stencil);
_stencil_val(u.x,0,0,o_stencil);
      
#line 299
_stencil_val_r(uf.z,0,0,0);  
    }


;; {        
       _stencil_val(u.z,0,-1,o_stencil);_stencil_val(u.z,0,0, o_stencil);_stencil_val(u.z,0,0, o_stencil); _stencil_val(u.z,0,1,o_stencil);_stencil_val(u.y,0,0, o_stencil);
_stencil_val(u.y,0,0,o_stencil);
      
#line 305
_stencil_val_r(uf.z,0,0,0);  
    }

;        

      


      

    
#line 308
_stencil_val_r(uf.z,0,0,0); 
  }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 292
foreach_face_generic(){is_face_x(){ {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (_const_fm.y && _const_fm.y) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }


    if (_const_fm.z && _const_fm.z) {
      double fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);
    }

    val(uf.x,0,0,0) *= _const_fm.x;
  }}end_is_face_x()
#line 292
is_face_y(){ {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (_const_fm.z && _const_fm.z) {
      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);
    }


    if (_const_fm.x && _const_fm.x) {
      double fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);
    }

    val(uf.y,0,0,0) *= _const_fm.y;
  }}end_is_face_y()
#line 292
is_face_z(){ {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;

    if (_const_fm.x && _const_fm.x) {
      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);
    }


    if (_const_fm.y && _const_fm.y) {
      double fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);
    }

    val(uf.z,0,0,0) *= _const_fm.z;
  }}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }

  delete ((scalar *)((vector[]){du,{{-1},{-1},{-1}}}));
}
#line 323
static int advection_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 323 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
      static int advection_term(const int i,const double t,Event *_ev){tracing("advection_term","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",323);
{
  if (!stokes) {
    prediction();
    mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);
    advection ((scalar *)((vector[]){u,{{-1},{-1},{-1}}}), uf, dt, (scalar *)((vector[]){g,{{-1},{-1},{-1}}}));
  }
}{end_tracing("advection_term","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",330);return 0;}end_tracing("advection_term","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",330);}







static void correction (double dt)
{
  foreach_stencil(1,{0},)
    {
      {_stencil_val(g.x,0,0,0);_stencil_val_r(u.x,0,0,0);  }
      
#line 342
{_stencil_val(g.y,0,0,0);_stencil_val_r(u.y,0,0,0);  }
      
#line 342
{_stencil_val(g.z,0,0,0);_stencil_val_r(u.z,0,0,0);  }}end_foreach_stencil()
   BEGIN_FOREACH{
#line 340
foreach()
    {
      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
      
#line 342
val(u.y,0,0,0) += dt*val(g.y,0,0,0);
      
#line 342
val(u.z,0,0,0) += dt*val(g.z,0,0,0);}end_foreach();}END_FOREACH 
}








static int viscous_term_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 352 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
      static int viscous_term(const int i,const double t,Event *_ev){tracing("viscous_term","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",352);
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax
#line 306 "/home/esteban/ProgramFile/basilisk/src/viscosity.h"
, NULL
#line 356 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
);
    correction (-dt);
  }




  if (!is_constant(a.x)) {
    vector af = a;
    trash (((vector[]){af,{{-1},{-1},{-1}}}));
    foreach_face_stencil(1,{0},){_stencil_is_face_x(){
      {_stencil_val_a(af.x,0,0,0);  }}end__stencil_is_face_x()
#line 366
_stencil_is_face_y(){
      {_stencil_val_a(af.y,0,0,0);  }}end__stencil_is_face_y()
#line 366
_stencil_is_face_z(){
      {_stencil_val_a(af.z,0,0,0);  }}end__stencil_is_face_z()}end_foreach_face_stencil()
     BEGIN_FOREACH{
#line 366
foreach_face_generic(){is_face_x(){
      val(af.x,0,0,0) = 0.;}end_is_face_x()
#line 366
is_face_y(){
      val(af.y,0,0,0) = 0.;}end_is_face_y()
#line 366
is_face_z(){
      val(af.z,0,0,0) = 0.;}end_is_face_z()}end_foreach_face_generic();}END_FOREACH 
  }
}{end_tracing("viscous_term","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",369);return 0;}end_tracing("viscous_term","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",369);}
#line 388
static int acceleration_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 388 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
      static int acceleration(const int i,const double t,Event *_ev){tracing("acceleration","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",388);
{
  trash (((vector[]){uf,{{-1},{-1},{-1}}}));
  if(!is_constant(fm.x) && !is_constant(a.x)){
  
#line 391
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val(a.x,0,0,0);_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 391
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val(a.y,0,0,0);_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()
#line 391
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);_stencil_val(a.z,0,0,0);_stencil_val_a(uf.z,0,0,0);    }}end__stencil_is_face_z()}end_foreach_face_stencil() BEGIN_FOREACH{
#line 391
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 391
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()
#line 391
is_face_z(){
    val(uf.z,0,0,0) = val(fm.z,0,0,0)*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*val(a.z,0,0,0));}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(is_constant(fm.x) && !is_constant(a.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 391
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);_stencil_val(a.x,0,0,0);_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 391
_stencil_is_face_y(){
    {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);_stencil_val(a.y,0,0,0);_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()
#line 391
_stencil_is_face_z(){
    {;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);_stencil_val(a.z,0,0,0);_stencil_val_a(uf.z,0,0,0);    }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 391
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val(a.x,0,0,0));}end_is_face_x()
#line 391
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val(a.y,0,0,0));}end_is_face_y()
#line 391
is_face_z(){
    val(uf.z,0,0,0) = _const_fm.z*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*val(a.z,0,0,0));}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(!is_constant(fm.x) && is_constant(a.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 391
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);;_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 391
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);;_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()
#line 391
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);;_stencil_val_a(uf.z,0,0,0);    }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 391
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = val(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()
#line 391
is_face_y(){
    val(uf.y,0,0,0) = val(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()
#line 391
is_face_z(){
    val(uf.z,0,0,0) = val(fm.z,0,0,0)*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*_const_a.z);}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 391
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {;_stencil_val(u.x,0,0,0); _stencil_val(u.x,0 -1,0,0);;_stencil_val_a(uf.x,0,0,0);    }}end__stencil_is_face_x()
#line 391
_stencil_is_face_y(){
    {;_stencil_val(u.y,0,0,0); _stencil_val(u.y,0,0 -1,0);;_stencil_val_a(uf.y,0,0,0);    }}end__stencil_is_face_y()
#line 391
_stencil_is_face_z(){
    {;_stencil_val(u.z,0,0,0); _stencil_val(u.z,0,0,0 -1);;_stencil_val_a(uf.z,0,0,0);    }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 391
foreach_face_generic(){is_face_x(){
    val(uf.x,0,0,0) = _const_fm.x*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*_const_a.x);}end_is_face_x()
#line 391
is_face_y(){
    val(uf.y,0,0,0) = _const_fm.y*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*_const_a.y);}end_is_face_y()
#line 391
is_face_z(){
    val(uf.z,0,0,0) = _const_fm.z*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*_const_a.z);}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }
}{end_tracing("acceleration","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",393);return 0;}end_tracing("acceleration","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",393);}
#line 402 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
void centered_gradient (scalar p, vector g)
{





  vector  gf=new_face_vector("gf");
  if(!is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){
  
#line 410
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(a.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 410
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(a.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 410
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);_stencil_val(a.z,0,0,0); _stencil_val(alpha.z,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil() BEGIN_FOREACH{
#line 410
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 410
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 410
is_face_z(){
    val(gf.z,0,0,0) = val(fm.z,0,0,0)*val(a.z,0,0,0) - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 410
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {;_stencil_val(a.x,0,0,0); _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 410
_stencil_is_face_y(){
    {;_stencil_val(a.y,0,0,0); _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 410
_stencil_is_face_z(){
    {;_stencil_val(a.z,0,0,0); _stencil_val(alpha.z,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 410
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 410
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 410
is_face_z(){
    val(gf.z,0,0,0) = _const_fm.z*val(a.z,0,0,0) - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(!is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 410
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);; _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 410
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);; _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 410
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);; _stencil_val(alpha.z,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 410
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 410
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 410
is_face_z(){
    val(gf.z,0,0,0) = val(fm.z,0,0,0)*_const_a.z - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);
  
#line 410
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {;; _stencil_val(alpha.x,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 410
_stencil_is_face_y(){
    {;; _stencil_val(alpha.y,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 410
_stencil_is_face_z(){
    {;; _stencil_val(alpha.z,0,0,0);_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 410
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - val(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 410
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - val(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 410
is_face_z(){
    val(gf.z,0,0,0) = _const_fm.z*_const_a.z - val(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(!is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 410
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);_stencil_val(a.x,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 410
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);_stencil_val(a.y,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 410
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);_stencil_val(a.z,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 410
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 410
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 410
is_face_z(){
    val(gf.z,0,0,0) = val(fm.z,0,0,0)*val(a.z,0,0,0) - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 410
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {;_stencil_val(a.x,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 410
_stencil_is_face_y(){
    {;_stencil_val(a.y,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 410
_stencil_is_face_z(){
    {;_stencil_val(a.z,0,0,0);;_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 410
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*val(a.x,0,0,0) - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 410
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*val(a.y,0,0,0) - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 410
is_face_z(){
    val(gf.z,0,0,0) = _const_fm.z*val(a.z,0,0,0) - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(!is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)){_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 410
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    { _stencil_val(fm.x,0,0,0);;;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 410
_stencil_is_face_y(){
    { _stencil_val(fm.y,0,0,0);;;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 410
_stencil_is_face_z(){
    { _stencil_val(fm.z,0,0,0);;;_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 410
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = val(fm.x,0,0,0)*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 410
is_face_y(){
    val(gf.y,0,0,0) = val(fm.y,0,0,0)*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 410
is_face_z(){
    val(gf.z,0,0,0) = val(fm.z,0,0,0)*_const_a.z - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 410
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {;;;_stencil_val(p,0,0,0); _stencil_val(p,-1,0,0);_stencil_val_a(gf.x,0,0,0);   }}end__stencil_is_face_x()
#line 410
_stencil_is_face_y(){
    {;;;_stencil_val(p,0,0,0); _stencil_val(p,0,-1,0);_stencil_val_a(gf.y,0,0,0);   }}end__stencil_is_face_y()
#line 410
_stencil_is_face_z(){
    {;;;_stencil_val(p,0,0,0); _stencil_val(p,0,0,-1);_stencil_val_a(gf.z,0,0,0);   }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 410
foreach_face_generic(){is_face_x(){
    val(gf.x,0,0,0) = _const_fm.x*_const_a.x - _const_alpha.x*(val(p,0,0,0) - val(p,-1,0,0))/Delta;}end_is_face_x()
#line 410
is_face_y(){
    val(gf.y,0,0,0) = _const_fm.y*_const_a.y - _const_alpha.y*(val(p,0,0,0) - val(p,0,-1,0))/Delta;}end_is_face_y()
#line 410
is_face_z(){
    val(gf.z,0,0,0) = _const_fm.z*_const_a.z - _const_alpha.z*(val(p,0,0,0) - val(p,0,0,-1))/Delta;}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }





  trash (((vector[]){g,{{-1},{-1},{-1}}}));
  if(!is_constant(fm.x)){
  
#line 418
foreach_stencil(1,{0},)
    {
      {_stencil_val(gf.x,0,0,0); _stencil_val(gf.x,1,0,0);_stencil_val(fm.x,0,0,0); _stencil_val(fm.x,1,0,0);_stencil_val_a(g.x,0,0,0);      }
      
#line 420
{_stencil_val(gf.y,0,0,0); _stencil_val(gf.y,0,1,0);_stencil_val(fm.y,0,0,0); _stencil_val(fm.y,0,1,0);_stencil_val_a(g.y,0,0,0);      }
      
#line 420
{_stencil_val(gf.z,0,0,0); _stencil_val(gf.z,0,0,1);_stencil_val(fm.z,0,0,0); _stencil_val(fm.z,0,0,1);_stencil_val_a(g.z,0,0,0);      }}end_foreach_stencil() BEGIN_FOREACH{
#line 418
foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val(fm.x,0,0,0) + val(fm.x,1,0,0) + 0.);
      
#line 420
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val(fm.y,0,0,0) + val(fm.y,0,1,0) + 0.);
      
#line 420
val(g.z,0,0,0) = (val(gf.z,0,0,0) + val(gf.z,0,0,1))/(val(fm.z,0,0,0) + val(fm.z,0,0,1) + 0.);}end_foreach();}END_FOREACH }else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 418
foreach_stencil(1,{0},)
    {
      {_stencil_val(gf.x,0,0,0); _stencil_val(gf.x,1,0,0);;;_stencil_val_a(g.x,0,0,0);      }
      
#line 420
{_stencil_val(gf.y,0,0,0); _stencil_val(gf.y,0,1,0);;;_stencil_val_a(g.y,0,0,0);      }
      
#line 420
{_stencil_val(gf.z,0,0,0); _stencil_val(gf.z,0,0,1);;;_stencil_val_a(g.z,0,0,0);      }}end_foreach_stencil()
   BEGIN_FOREACH{
#line 418
foreach()
    {
      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(_const_fm.x + _const_fm.x + 0.);
      
#line 420
val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(_const_fm.y + _const_fm.y + 0.);
      
#line 420
val(g.z,0,0,0) = (val(gf.z,0,0,0) + val(gf.z,0,0,1))/(_const_fm.z + _const_fm.z + 0.);}end_foreach();}END_FOREACH }delete((scalar*)((vector[]){gf,{{-1},{-1},{-1}}}));
}






static int projection_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}







#line 428
      static int projection(const int i,const double t,Event *_ev){tracing("projection","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",428);
{
  mgp = project (uf, p, alpha, dt, mgp.nrelax);
  centered_gradient (p, g);




  correction (dt);
}{end_tracing("projection","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",437);return 0;}end_tracing("projection","/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",437);}





static int end_timestep_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}






#line 443
static int end_timestep(const int i,const double t,Event *_ev){;return 0;}
#line 12 "main3D.c"

#line 1 "two-phase.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/two-phase.h"
#line 13 "/home/esteban/ProgramFile/basilisk/src/two-phase.h"
#line 1 "vof.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/vof.h"
#line 27 "/home/esteban/ProgramFile/basilisk/src/vof.h"





#line 1 "fractions.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
#line 12 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/geometry.h"
#line 35 "/home/esteban/ProgramFile/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    do { double _tmp_ = n1; n1 = n2; n2 = _tmp_; } while(false);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}



double plane_alpha (double c, coord n)
{
  double alpha;
  coord n1;

  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);

  double m1, m2, m3;
  m1 = min(n1.x, n1.y);
  m3 = max(n1.x, n1.y);
  m2 = n1.z;
  if (m2 < m1) {
    double tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    double tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  double m12 = m1 + m2;
  double pr = max(6.*m1*m2*m3, 1e-50);
  double V1 = m1*m1*m1/pr;
  double V2 = V1 + (m2 - m1)/(2.*m3), V3;
  double mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  c = clamp (c, 0., 1.);
  double ch = min(c, 1. - c);
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    double p12 = sqrt (2.*m1*m2);
    double q = 3.*(m12 - 2.*m3*ch)/(4.*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 <= m3)
    alpha = m3*ch + mm/2.;
  else {
    double p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);
    double q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;
  if (n.z < 0.)
    alpha += n.z;

  return alpha - (n.x + n.y + n.z)/2.;;
}
#line 163 "/home/esteban/ProgramFile/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}



double plane_volume (coord n, double alpha)
{
  double al = alpha + (n.x + n.y + n.z)/2. +
    max(0., -n.x) + max(0., -n.y) + max(0., -n.z);
  if (al <= 0.)
    return 0.;
  double tmp = fabs(n.x) + fabs(n.y) + fabs(n.z);
  if (al >= tmp)
    return 1.;
  if (tmp < 1e-10)
    return 0.;
  double n1 = fabs(n.x)/tmp;
  double n2 = fabs(n.y)/tmp;
  double n3 = fabs(n.z)/tmp;
  al = max(0., min(1., al/tmp));
  double al0 = min(al, 1. - al);
  double b1 = min(n1, n2);
  double b3 = max(n1, n2);
  double b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  double b12 = b1 + b2;
  double bm = min(b12, b3);
  double pr = max(6.*b1*b2*b3, 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) + b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  double volume = al <= 0.5 ? tmp : 1. - tmp;
  return clamp (volume, 0., 1.);
}
#line 267 "/home/esteban/ProgramFile/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
   {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  } 
#line 270
{
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  } 
#line 270
{
    alpha -= n.z*(b.z + a.z)/2.;
    n1.z = n.z*(b.z - a.z);
  }
  return plane_volume (n1, alpha);
}
#line 307 "/home/esteban/ProgramFile/basilisk/src/geometry.h"
static coord cube_edge[12][2] = {
  {{0.,0.,0.},{1.,0.,0.}},{{0.,0.,1.},{1.,0.,1.}},
  {{0.,1.,1.},{1.,1.,1.}},{{0.,1.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,1.,0.}},{{0.,0.,1.},{0.,1.,1.}},
  {{1.,0.,1.},{1.,1.,1.}},{{1.,0.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,0.,1.}},{{1.,0.,0.},{1.,0.,1.}},
  {{1.,1.,0.},{1.,1.,1.}},{{0.,1.,0.},{0.,1.,1.}}
};




static int cube_connect[12][2][4] = {
  {{9, 1, 8}, {4, 3, 7}},
  {{6, 2, 5}, {8, 0, 9}},
  {{10, 3, 11}, {5, 1, 6}},
  {{7, 0, 4}, {11, 2, 10}},
  {{3, 7, 0}, {8, 5, 11}},
  {{11, 4, 8}, {1, 6, 2}},
  {{2, 5, 1}, {9, 7, 10}},
  {{10, 6, 9}, {0, 4, 3}},
  {{5, 11, 4}, {0, 9, 1}},
  {{1, 8, 0}, {7, 10, 6}},
  {{6, 9, 7}, {3, 11, 2}},
  {{2, 10, 3}, {4, 8, 5}}
};

int facets (coord n, double alpha, coord v[12], double h)
{
  coord a[12];
  int orient[12];

  for (int i = 0; i < 12; i++) {
    coord e, d;
    double den = 0., t = alpha;
     {
      d.x = h*(cube_edge[i][0].x - 0.5);
      e.x = h*(cube_edge[i][1].x - 0.5);
      den += n.x*(e.x - d.x);
      t -= n.x*d.x;
    } 
#line 342
{
      d.y = h*(cube_edge[i][0].y - 0.5);
      e.y = h*(cube_edge[i][1].y - 0.5);
      den += n.y*(e.y - d.y);
      t -= n.y*d.y;
    } 
#line 342
{
      d.z = h*(cube_edge[i][0].z - 0.5);
      e.z = h*(cube_edge[i][1].z - 0.5);
      den += n.z*(e.z - d.z);
      t -= n.z*d.z;
    }
    orient[i] = -1;
    if (fabs (den) > 1e-10) {
      t /= den;
      if (t >= 0. && t < 1.) {
 double s = - alpha;
  {
   a[i].x = d.x + t*(e.x - d.x);
   s += n.x*e.x;
 } 
#line 353
{
   a[i].y = d.y + t*(e.y - d.y);
   s += n.y*e.y;
 } 
#line 353
{
   a[i].z = d.z + t*(e.z - d.z);
   s += n.z*e.z;
 }
 orient[i] = (s > 0.);
      }
    }
  }

  for (int i = 0; i < 12; i++) {
    int nv = 0, e = i;
    while (orient[e] >= 0) {
      int m = 0, * ne = cube_connect[e][orient[e]];
      v[nv++] = a[e];
      orient[e] = -1;
      while (m < 3 && orient[e] < 0)
 e = ne[m++];
    }
    if (nv > 2)
      return nv;
  }
  return 0;
}






double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 388
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->x = p->y = p->z = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
    
#line 399
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

   {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  } 
#line 424
{
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }

  return sqrt (ax*ax + ay*ay);
}




double plane_area_center (coord m, double alpha, coord * p)
{
  
    if (fabs (m.x) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.y;
      ((double *)&n)[1] = m.z;
      double length = line_length_center (n, alpha, &q);
      p->x = 0.;
      p->y = ((double *)&q)[0];
      p->z = ((double *)&q)[1];
      return length;
    }
    
#line 441
if (fabs (m.y) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.z;
      ((double *)&n)[1] = m.x;
      double length = line_length_center (n, alpha, &q);
      p->y = 0.;
      p->z = ((double *)&q)[0];
      p->x = ((double *)&q)[1];
      return length;
    }
    
#line 441
if (fabs (m.z) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.x;
      ((double *)&n)[1] = m.y;
      double length = line_length_center (n, alpha, &q);
      p->z = 0.;
      p->x = ((double *)&q)[0];
      p->y = ((double *)&q)[1];
      return length;
    }

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 455
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }
    
#line 455
if (n.z < 0.) {
      alpha -= n.z;
      n.z = - n.z;
    }

  double amax = n.x + n.y + n.z;
  if (alpha < 0. || alpha > amax) {
    p->x = p->y = p->z = 0.;
    return 0.;
  }

  double area = sq(alpha);
  p->x = p->y = p->z = area*alpha;

   {
    double b = alpha - n.x;
    if (b > 0.) {
      area -= b*b;
      p->x -= b*b*(2.*n.x + alpha);
      p->y -= b*b*b;
      p->z -= b*b*b;
    }
  } 
#line 469
{
    double b = alpha - n.y;
    if (b > 0.) {
      area -= b*b;
      p->y -= b*b*(2.*n.y + alpha);
      p->z -= b*b*b;
      p->x -= b*b*b;
    }
  } 
#line 469
{
    double b = alpha - n.z;
    if (b > 0.) {
      area -= b*b;
      p->z -= b*b*(2.*n.z + alpha);
      p->x -= b*b*b;
      p->y -= b*b*b;
    }
  }

  amax = alpha - amax;
   {
    double b = amax + n.x;
    if (b > 0.) {
      area += b*b;
      p->y += b*b*(2.*n.y + alpha - n.z);
      p->z += b*b*(2.*n.z + alpha - n.y);
      p->x += b*b*b;
    }
  } 
#line 480
{
    double b = amax + n.y;
    if (b > 0.) {
      area += b*b;
      p->z += b*b*(2.*n.z + alpha - n.x);
      p->x += b*b*(2.*n.x + alpha - n.z);
      p->y += b*b*b;
    }
  } 
#line 480
{
    double b = amax + n.z;
    if (b > 0.) {
      area += b*b;
      p->x += b*b*(2.*n.x + alpha - n.y);
      p->y += b*b*(2.*n.y + alpha - n.x);
      p->z += b*b*b;
    }
  }

  area *= 3.;
   {
    if (area) {
      p->x /= area*n.x;
      p->x = clamp (p->x, 0., 1.);
    }
    else
      p->x = 0.;
    if (m.x < 0.) p->x = 1. - p->x;
    p->x -= 0.5;
  } 
#line 491
{
    if (area) {
      p->y /= area*n.y;
      p->y = clamp (p->y, 0., 1.);
    }
    else
      p->y = 0.;
    if (m.y < 0.) p->y = 1. - p->y;
    p->y -= 0.5;
  } 
#line 491
{
    if (area) {
      p->z /= area*n.z;
      p->z = clamp (p->z, 0., 1.);
    }
    else
      p->z = 0.;
    if (m.z < 0.) p->z = 1. - p->z;
    p->z -= 0.5;
  }

  return area*sqrt (1./(sq(n.x)*sq(n.y)) +
      1./(sq(n.x)*sq(n.z)) +
      1./(sq(n.z)*sq(n.y)))/6.;
}






void line_center (coord m, double alpha, double a, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 518
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = -0.5;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return;
  }

  
    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = sign(m.y)*(a/2. - 0.5);
      return;
    }
    
#line 535
if (n.y < 1e-4) {
      p->y = 0.;
      p->x = sign(m.x)*(a/2. - 0.5);
      return;
    }

  p->x = p->y = cube(alpha);

   {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= sq(b)*(alpha + 2.*n.x);
      p->y -= cube(b);
    }
  } 
#line 543
{
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= sq(b)*(alpha + 2.*n.y);
      p->x -= cube(b);
    }
  }

   {
    p->x /= 6.*sq(n.x)*n.y*a;
    p->x = sign(m.x)*(p->x - 0.5);
  } 
#line 551
{
    p->y /= 6.*sq(n.y)*n.x*a;
    p->y = sign(m.y)*(p->y - 0.5);
  }
}
#line 564 "/home/esteban/ProgramFile/basilisk/src/geometry.h"
void plane_center (coord m, double alpha, double a, coord * p)
{
  
    if (fabs (m.x) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.y;
      ((double *)&n)[1] = m.z;
      line_center (n, alpha, a, &q);
      p->x = 0.;
      p->y = ((double *)&q)[0];
      p->z = ((double *)&q)[1];
      return;
    }
    
#line 567
if (fabs (m.y) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.z;
      ((double *)&n)[1] = m.x;
      line_center (n, alpha, a, &q);
      p->y = 0.;
      p->z = ((double *)&q)[0];
      p->x = ((double *)&q)[1];
      return;
    }
    
#line 567
if (fabs (m.z) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.x;
      ((double *)&n)[1] = m.y;
      line_center (n, alpha, a, &q);
      p->z = 0.;
      p->x = ((double *)&q)[0];
      p->y = ((double *)&q)[1];
      return;
    }

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
    
#line 581
if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }
    
#line 581
if (n.z < 0.) {
      alpha -= n.z;
      n.z = - n.z;
    }

  if (alpha <= 0. || a == 0.) {
    p->x = p->y = p->z = -0.5;
    return;
  }

  if (alpha >= n.x + n.y + n.z || a == 1.) {
    p->x = p->y = p->z = 0.;
    return;
  }

  p->x = p->y = p->z = sq(sq(alpha));
   {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= cube(b)*(3.*n.x + alpha);
      p->y -= sq(sq(b));
      p->z -= sq(sq(b));
    }
  } 
#line 597
{
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= cube(b)*(3.*n.y + alpha);
      p->z -= sq(sq(b));
      p->x -= sq(sq(b));
    }
  } 
#line 597
{
    double b = alpha - n.z;
    if (b > 0.) {
      p->z -= cube(b)*(3.*n.z + alpha);
      p->x -= sq(sq(b));
      p->y -= sq(sq(b));
    }
  }

  double amax = alpha - (n.x + n.y + n.z);
   {
    double b = amax + n.z;
    if (b > 0.) {
      p->x += cube(b)*(3.*n.x + alpha - n.y);
      p->y += cube(b)*(3.*n.y + alpha - n.x);
      p->z += sq(sq(b));
    }
  } 
#line 607
{
    double b = amax + n.x;
    if (b > 0.) {
      p->y += cube(b)*(3.*n.y + alpha - n.z);
      p->z += cube(b)*(3.*n.z + alpha - n.y);
      p->x += sq(sq(b));
    }
  } 
#line 607
{
    double b = amax + n.y;
    if (b > 0.) {
      p->z += cube(b)*(3.*n.z + alpha - n.x);
      p->x += cube(b)*(3.*n.x + alpha - n.z);
      p->y += sq(sq(b));
    }
  }

  double b = 24.*n.x*n.y*n.z*a;
   {
    p->x /= b*n.x;
    p->x = sign(m.x)*(p->x - 0.5);
  } 
#line 617
{
    p->y /= b*n.y;
    p->y = sign(m.y)*(p->y - 0.5);
  } 
#line 617
{
    p->z /= b*n.z;
    p->z = sign(m.z)*(p->z - 0.5);
  }
}
#line 13 "/home/esteban/ProgramFile/basilisk/src/fractions.h"







#line 1 "myc.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/myc.h"
#line 16 "/home/esteban/ProgramFile/basilisk/src/myc.h"
coord mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double m1,m2,m[4][3],t0,t1,t2;
  int cn;



  m1 = val(c,-1,0,-1) + val(c,-1,0,1) + val(c,-1,-1,0) + val(c,-1,1,0) +
       val(c,-1,0,0);
  m2 = val(c,1,0,-1) + val(c,1,0,1) + val(c,1,-1,0) + val(c,1,1,0) +
       val(c,1,0,0);
  m[0][0] = m1 > m2 ? 1. : -1.;

  m1 = val(c,-1,-1,0)+ val(c,1,-1,0)+ val(c,0,-1,0);
  m2 = val(c,-1,1,0)+ val(c,1,1,0)+ val(c,0,1,0);
  m[0][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1)+ val(c,1,0,-1)+ val(c,0,0,-1);
  m2 = val(c,-1,0,1)+ val(c,1,0,1)+ val(c,0,0,1);
  m[0][2] = 0.5*(m1-m2);



  m1 = val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,0);
  m2 = val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,0);
  m[1][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1) + val(c,0,-1,1) + val(c,1,-1,0) + val(c,-1,-1,0) +
       val(c,0,-1,0);
  m2 = val(c,0,1,-1) + val(c,0,1,1) + val(c,1,1,0) + val(c,-1,1,0) +
       val(c,0,1,0);
  m[1][1] = m1 > m2 ? 1. : -1.;

  m1 = val(c,0,-1,-1)+ val(c,0,0,-1)+ val(c,0,1,-1);
  m2 = val(c,0,-1,1)+ val(c,0,0,1)+ val(c,0,1,1);
  m[1][2] = 0.5*(m1-m2);




  m1 = val(c,-1,0,-1)+ val(c,-1,0,1)+ val(c,-1,0,0);
  m2 = val(c,1,0,-1)+ val(c,1,0,1)+ val(c,1,0,0);
  m[2][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1)+ val(c,0,-1,1)+ val(c,0,-1,0);
  m2 = val(c,0,1,-1)+ val(c,0,1,1)+ val(c,0,1,0);
  m[2][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1) +
       val(c,0,0,-1);
  m2 = val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1) +
       val(c,0,0,1);
  m[2][2] = m1 > m2 ? 1. : -1.;


  t0 = fabs(m[0][0]) + fabs(m[0][1]) + fabs(m[0][2]);
  m[0][0] /= t0;
  m[0][1] /= t0;
  m[0][2] /= t0;

  t0 = fabs(m[1][0]) + fabs(m[1][1]) + fabs(m[1][2]);
  m[1][0] /= t0;
  m[1][1] /= t0;
  m[1][2] /= t0;

  t0 = fabs(m[2][0]) + fabs(m[2][1]) + fabs(m[2][2]);
  m[2][0] /= t0;
  m[2][1] /= t0;
  m[2][2] /= t0;


  t0 = fabs(m[0][0]);
  t1 = fabs(m[1][1]);
  t2 = fabs(m[2][2]);

  cn = 0;
  if (t1 > t0) {
    t0 = t1;
    cn = 1;
  }
  if (t2 > t0)
    cn = 2;


  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,-1,-1,1) + val(c,-1,1,1) +
       2.*(val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,-1) + val(c,-1,0,1)) +
       4.*val(c,-1,0,0);
  m2 = val(c,1,-1,-1) + val(c,1,1,-1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,-1) + val(c,1,0,1)) +
       4.*val(c,1,0,0);
  m[3][0] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,-1,1) + val(c,1,-1,-1) + val(c,1,-1,1) +
       2.*( val(c,-1,-1,0) + val(c,1,-1,0) + val(c,0,-1,-1) + val(c,0,-1,1)) +
       4.*val(c,0,-1,0);
  m2 = val(c,-1,1,-1) + val(c,-1,1,1) + val(c,1,1,-1) + val(c,1,1,1) +
       2.*(val(c,-1,1,0) + val(c,1,1,0) + val(c,0,1,-1) + val(c,0,1,1)) +
       4.*val(c,0,1,0);
  m[3][1] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,1,-1,-1) + val(c,1,1,-1) +
       2.*(val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1)) +
       4.*val(c,0,0,-1);
  m2 = val(c,-1,-1,1) + val(c,-1,1,1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1)) +
       4.*val(c,0,0,1);
  m[3][2] = m1 - m2;


  t0 = fabs(m[3][0]) + fabs(m[3][1]) + fabs(m[3][2]);
  if (t0 < 1e-30) {
    coord mxyz = {1., 0., 0.};
    return mxyz;
  }

  m[3][0] /= t0;
  m[3][1] /= t0;
  m[3][2] /= t0;


  t0 = fabs (m[3][0]);
  t1 = fabs (m[3][1]);
  t2 = fabs (m[3][2]);
  if (t1 > t0)
    t0 = t1;
  if (t2 > t0)
    t0 = t2;

  if (fabs(m[cn][cn]) > t0)
    cn = 3;


  coord mxyz = {m[cn][0], m[cn][1], m[cn][2]};
  return mxyz;
}
#line 13 "/home/esteban/ProgramFile/basilisk/src/fractions.h"







#line 1 "myc.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/myc.h"
#line 16 "/home/esteban/ProgramFile/basilisk/src/myc.h"
static void _stencil_mycs (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;  
  
   



_stencil_val(c,-1,0,-1); _stencil_val(c,-1,0,1); _stencil_val(c,-1,-1,0); _stencil_val(c,-1,1,0);
       _stencil_val(c,-1,0,0); 



       
#line 25
_stencil_val(c,1,0,-1); _stencil_val(c,1,0,1); _stencil_val(c,1,-1,0); _stencil_val(c,1,1,0);
       _stencil_val(c,1,0,0);  
      
          
_stencil_val(c,-1,-1,0); _stencil_val(c,1,-1,0); _stencil_val(c,0,-1,0);

    
#line 30
_stencil_val(c,-1,1,0); _stencil_val(c,1,1,0); _stencil_val(c,0,1,0);
   
     
_stencil_val(c,-1,0,-1); _stencil_val(c,1,0,-1); _stencil_val(c,0,0,-1);

    
#line 34
_stencil_val(c,-1,0,1); _stencil_val(c,1,0,1); _stencil_val(c,0,0,1);
   
     


_stencil_val(c,-1,-1,0); _stencil_val(c,-1,1,0); _stencil_val(c,-1,0,0);



      
#line 40
_stencil_val(c,1,-1,0); _stencil_val(c,1,1,0); _stencil_val(c,1,0,0);
     
     
_stencil_val(c,0,-1,-1); _stencil_val(c,0,-1,1); _stencil_val(c,1,-1,0); _stencil_val(c,-1,-1,0);
       _stencil_val(c,0,-1,0);

        
#line 45
_stencil_val(c,0,1,-1); _stencil_val(c,0,1,1); _stencil_val(c,1,1,0); _stencil_val(c,-1,1,0);
       _stencil_val(c,0,1,0);
       
           
_stencil_val(c,0,-1,-1); _stencil_val(c,0,0,-1); _stencil_val(c,0,1,-1);

    
#line 50
_stencil_val(c,0,-1,1); _stencil_val(c,0,0,1); _stencil_val(c,0,1,1);
   
     



_stencil_val(c,-1,0,-1); _stencil_val(c,-1,0,1); _stencil_val(c,-1,0,0);




    
#line 57
_stencil_val(c,1,0,-1); _stencil_val(c,1,0,1); _stencil_val(c,1,0,0);
   
     
_stencil_val(c,0,-1,-1); _stencil_val(c,0,-1,1); _stencil_val(c,0,-1,0);

    
#line 61
_stencil_val(c,0,1,-1); _stencil_val(c,0,1,1); _stencil_val(c,0,1,0);
   
     
_stencil_val(c,-1,0,-1); _stencil_val(c,1,0,-1); _stencil_val(c,0,-1,-1); _stencil_val(c,0,1,-1);
       _stencil_val(c,0,0,-1);

        
#line 66
_stencil_val(c,-1,0,1); _stencil_val(c,1,0,1); _stencil_val(c,0,-1,1); _stencil_val(c,0,1,1);
       _stencil_val(c,0,0,1);    
       
          


       
    
    
    

        
    
    
    

        
    
    
    


    
   
   

    
      
     
   
      
     
      

_stencil_val(c,-1,-1,-1); _stencil_val(c,-1,1,-1); _stencil_val(c,-1,-1,1); _stencil_val(c,-1,1,1);
_stencil_val(c,-1,-1,0); _stencil_val(c,-1,1,0); _stencil_val(c,-1,0,-1); _stencil_val(c,-1,0,1);
_stencil_val(c,-1,0,0);


       
           
        
#line 103
_stencil_val(c,1,-1,-1); _stencil_val(c,1,1,-1); _stencil_val(c,1,-1,1); _stencil_val(c,1,1,1);
_stencil_val(c,1,-1,0); _stencil_val(c,1,1,0); _stencil_val(c,1,0,-1); _stencil_val(c,1,0,1);
_stencil_val(c,1,0,0);
       
           
       
       
#line 108
_stencil_val(c,-1,-1,-1); _stencil_val(c,-1,-1,1); _stencil_val(c,1,-1,-1); _stencil_val(c,1,-1,1); 
_stencil_val(c,-1,-1,0); _stencil_val(c,1,-1,0); _stencil_val(c,0,-1,-1); _stencil_val(c,0,-1,1);
_stencil_val(c,0,-1,0);

       
           
        
#line 111
_stencil_val(c,-1,1,-1); _stencil_val(c,-1,1,1); _stencil_val(c,1,1,-1); _stencil_val(c,1,1,1);
_stencil_val(c,-1,1,0); _stencil_val(c,1,1,0); _stencil_val(c,0,1,-1); _stencil_val(c,0,1,1);
_stencil_val(c,0,1,0);
       
           
       
       
#line 116
_stencil_val(c,-1,-1,-1); _stencil_val(c,-1,1,-1); _stencil_val(c,1,-1,-1); _stencil_val(c,1,1,-1);
_stencil_val(c,-1,0,-1); _stencil_val(c,1,0,-1); _stencil_val(c,0,-1,-1); _stencil_val(c,0,1,-1);
_stencil_val(c,0,0,-1);

       
           
        
#line 119
_stencil_val(c,-1,-1,1); _stencil_val(c,-1,1,1); _stencil_val(c,1,-1,1); _stencil_val(c,1,1,1);
_stencil_val(c,-1,0,1); _stencil_val(c,1,0,1); _stencil_val(c,0,-1,1); _stencil_val(c,0,1,1);
_stencil_val(c,0,0,1);          
     
    
          
       
           
       
      


        
     

    
    
    


     
     
     
  
      
  
      
     

     


  
  
#line 149
return ;
}
#line 21 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
#line 120 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
     
void fractions (scalar Phi, scalar c,
  vector s, double val)
{tracing("fractions","/home/esteban/ProgramFile/basilisk/src/fractions.h",121);

  vector   as=(s).x.i>0?(s):new_face_vector("as");
#line 134 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
  vector  p=new_vector("p");
#line 146 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
   BEGIN_FOREACH{foreach_vertex_aux() { if (_a.x < point.n + 2) {





    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {






      val(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < val)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 171 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);
  } 
#line 146
if (_a.y < point.n + 2) {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {






      val(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < val)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 171 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);
  } 
#line 146
if (_a.z < point.n + 2) {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,0,1) - val) < 0.) {






      val(p.z,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,0,1));
      if (val(Phi,0,0,0) < val)
 val(p.z,0,0,0) = 1. - val(p.z,0,0,0);
    }
#line 171 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
    else
      val(p.z,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,0,1) > val);
  }}end_foreach_vertex_aux();}END_FOREACH 
#line 190 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
  
    _attribute[p.x.i].dirty = false;
    
#line 191
_attribute[p.y.i].dirty = false;
    
#line 191
_attribute[p.z.i].dirty = false;

  scalar s_x = as.x, s_y = as.y, s_z = as.z;
  foreach_face_stencil(1,{0},){_stencil_is_face_z(){




  {    
#line 231 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.y,0,0,0); _stencil_val(p.y,1,0,0);  
       
       
    
#line 236
} 
#line 233
{ 
_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,1,0);  
       
       
    
#line 236
}





{
      { _stencil_val(p.x,0,0,0);_stencil_val_a(s_z,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.x,0,i,0); _stencil_val(p.x,0,i,0); {       
     _stencil_val(p.x,0,i,0);_stencil_val(Phi,0,i,0); 
          
     
   }      }
   
#line 261
{_stencil_val(p.y,i,0,0); _stencil_val(p.y,i,0,0); {       
     _stencil_val(p.y,i,0,0);_stencil_val(Phi,i,0,0); 
          
     
   }      }}








{
 {_stencil_val(p.x,0,0,0); _stencil_val(p.y,0,0,0);_stencil_val_a(s_z,0,0,0);   }
{
 {_stencil_val_a(s_z,0,0,0);     } 
{

_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,1,0); _stencil_val(p.y,0,0,0); _stencil_val(p.y,1,0,0);

 
#line 280
_stencil_val_a(s_z,0,0,0);       



      }}}
#line 274 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
         
          
      
    







}}





       
    
  
#line 286
}}end__stencil_is_face_z()
#line 194
_stencil_is_face_x(){




  {    
#line 231 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.z,0,0,0); _stencil_val(p.z,0,1,0);  
       
       
    
#line 236
} 
#line 233
{ 
_stencil_val(p.y,0,0,0); _stencil_val(p.y,0,0,1);  
       
       
    
#line 236
}





{
      { _stencil_val(p.y,0,0,0);_stencil_val_a(s_x,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.y,0,0,i); _stencil_val(p.y,0,0,i); {       
     _stencil_val(p.y,0,0,i);_stencil_val(Phi,0,0,i); 
          
     
   }      }
   
#line 261
{_stencil_val(p.z,0,i,0); _stencil_val(p.z,0,i,0); {       
     _stencil_val(p.z,0,i,0);_stencil_val(Phi,0,i,0); 
          
     
   }      }}








{
 {_stencil_val(p.y,0,0,0); _stencil_val(p.z,0,0,0);_stencil_val_a(s_x,0,0,0);   }
{
 {_stencil_val_a(s_x,0,0,0);     } 
{

_stencil_val(p.y,0,0,0); _stencil_val(p.y,0,0,1); _stencil_val(p.z,0,0,0); _stencil_val(p.z,0,1,0);

 
#line 280
_stencil_val_a(s_x,0,0,0);       



      }}}
#line 274 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
         
          
      
    







}}





       
    
  
#line 286
}}end__stencil_is_face_x()
#line 194
_stencil_is_face_y(){




  {    
#line 231 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
    
    
     { 
_stencil_val(p.x,0,0,0); _stencil_val(p.x,0,0,1);  
       
       
    
#line 236
} 
#line 233
{ 
_stencil_val(p.z,0,0,0); _stencil_val(p.z,1,0,0);  
       
       
    
#line 236
}





{
      { _stencil_val(p.z,0,0,0);_stencil_val_a(s_y,0,0,0); } 
{      





      
   






      
      for (int i = 0; i <= 1; i++)
 {
   {_stencil_val(p.z,i,0,0); _stencil_val(p.z,i,0,0); {       
     _stencil_val(p.z,i,0,0);_stencil_val(Phi,i,0,0); 
          
     
   }      }
   
#line 261
{_stencil_val(p.x,0,0,i); _stencil_val(p.x,0,0,i); {       
     _stencil_val(p.x,0,0,i);_stencil_val(Phi,0,0,i); 
          
     
   }      }}








{
 {_stencil_val(p.z,0,0,0); _stencil_val(p.x,0,0,0);_stencil_val_a(s_y,0,0,0);   }
{
 {_stencil_val_a(s_y,0,0,0);     } 
{

_stencil_val(p.z,0,0,0); _stencil_val(p.z,1,0,0); _stencil_val(p.x,0,0,0); _stencil_val(p.x,0,0,1);

 
#line 280
_stencil_val_a(s_y,0,0,0);       



      }}}
#line 274 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
         
          
      
    







}}





       
    
  
#line 286
}}end__stencil_is_face_y()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 194
foreach_face_generic(){is_face_z(){




  {
#line 231 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    } 
#line 233
{
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      
 n.x /= nn;
 
#line 251
n.y /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
   
#line 261
if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 274 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));
      else if (ni != 4)
 val(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
      else {

 val(s_z,0,0,0) = (val(p.x,0,0,0) + val(p.x,0,1,0) + val(p.y,0,0,0) + val(p.y,1,0,0) > 2.);



      }
    }
  }}end_is_face_z()
#line 194
is_face_x(){




  {
#line 231 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.y = val(p.z,0,0,0) - val(p.z,0,1,0);
      nn += fabs(n.y);
    } 
#line 233
{
      n.z = val(p.y,0,0,0) - val(p.y,0,0,1);
      nn += fabs(n.z);
    }





    if (nn == 0.)
      val(s_x,0,0,0) = val(p.y,0,0,0);
    else {





      
 n.y /= nn;
 
#line 251
n.z /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.y,0,0,i) > 0. && val(p.y,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i) - val)*(val(p.y,0,0,i) - 0.5);
     alpha += n.y*a + n.z*(i - 0.5);
     ni++;
   }
   
#line 261
if (val(p.z,0,i,0) > 0. && val(p.z,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.z,0,i,0) - 0.5);
     alpha += n.z*a + n.y*(i - 0.5);
     ni++;
   }}
#line 274 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_x,0,0,0) = max (val(p.y,0,0,0), val(p.z,0,0,0));
      else if (ni != 4)
 val(s_x,0,0,0) = line_area (n.y, n.z, alpha/ni);
      else {

 val(s_x,0,0,0) = (val(p.y,0,0,0) + val(p.y,0,0,1) + val(p.z,0,0,0) + val(p.z,0,1,0) > 2.);



      }
    }
  }}end_is_face_x()
#line 194
is_face_y(){




  {
#line 231 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
     {
      n.z = val(p.x,0,0,0) - val(p.x,0,0,1);
      nn += fabs(n.z);
    } 
#line 233
{
      n.x = val(p.z,0,0,0) - val(p.z,1,0,0);
      nn += fabs(n.x);
    }





    if (nn == 0.)
      val(s_y,0,0,0) = val(p.z,0,0,0);
    else {





      
 n.z /= nn;
 
#line 251
n.x /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
   if (val(p.z,i,0,0) > 0. && val(p.z,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.z,i,0,0) - 0.5);
     alpha += n.z*a + n.x*(i - 0.5);
     ni++;
   }
   
#line 261
if (val(p.x,0,0,i) > 0. && val(p.x,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i) - val)*(val(p.x,0,0,i) - 0.5);
     alpha += n.x*a + n.z*(i - 0.5);
     ni++;
   }}
#line 274 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_y,0,0,0) = max (val(p.z,0,0,0), val(p.x,0,0,0));
      else if (ni != 4)
 val(s_y,0,0,0) = line_area (n.z, n.x, alpha/ni);
      else {

 val(s_y,0,0,0) = (val(p.z,0,0,0) + val(p.z,1,0,0) + val(p.x,0,0,0) + val(p.x,0,0,1) > 2.);



      }
    }
  }}end_is_face_y()}end_foreach_face_generic();}END_FOREACH 







  foreach_stencil(1,{0},) {    




    
    
     { 
_stencil_val(as.x,0,0,0); _stencil_val(as.x,1,0,0);  
       
       
    
#line 304
} 
#line 301
{ 
_stencil_val(as.y,0,0,0); _stencil_val(as.y,0,1,0);  
       
       
    
#line 304
} 
#line 301
{ 
_stencil_val(as.z,0,0,0); _stencil_val(as.z,0,0,1);  
       
       
    
#line 304
}
{
      { _stencil_val(as.x,0,0,0);_stencil_val_a(c,0,0,0); } 
{      
      
   






      
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   {
     {_stencil_val(p.x,0,i,j); _stencil_val(p.x,0,i,j); {       
       _stencil_val(p.x,0,i,j);_stencil_val(Phi,0,i,j); 
                
       
     }      }
     
#line 320
{_stencil_val(p.y,j,0,i); _stencil_val(p.y,j,0,i); {       
       _stencil_val(p.y,j,0,i);_stencil_val(Phi,j,0,i); 
                
       
     }      }
     
#line 320
{_stencil_val(p.z,i,j,0); _stencil_val(p.z,i,j,0); {       
       _stencil_val(p.z,i,j,0);_stencil_val(Phi,i,j,0); 
                
       
     }      }}




{
 { _stencil_val(as.x,0,0,0);_stencil_val_a(c,0,0,0); }
{
 {_stencil_val_a(c,0,0,0);  }
 
{_stencil_val_a(c,0,0,0);    }}}




         
              
      
    
#line 335
}}
       
    
  
#line 336
}end_foreach_stencil()







   BEGIN_FOREACH{
#line 294
foreach() {




    coord n;
    double nn = 0.;
     {
      n.x = val(as.x,0,0,0) - val(as.x,1,0,0);
      nn += fabs(n.x);
    } 
#line 301
{
      n.y = val(as.y,0,0,0) - val(as.y,0,1,0);
      nn += fabs(n.y);
    } 
#line 301
{
      n.z = val(as.z,0,0,0) - val(as.z,0,0,1);
      nn += fabs(n.z);
    }
    if (nn == 0.)
      val(c,0,0,0) = val(as.x,0,0,0);
    else {
      
 n.x /= nn;
 
#line 309
n.y /= nn;
 
#line 309
n.z /= nn;






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   {
     if (val(p.x,0,i,j) > 0. && val(p.x,0,i,j) < 1.) {
       double a = sign(val(Phi,0,i,j) - val)*(val(p.x,0,i,j) - 0.5);
       alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);
       ni++;
     }
     
#line 320
if (val(p.y,j,0,i) > 0. && val(p.y,j,0,i) < 1.) {
       double a = sign(val(Phi,j,0,i) - val)*(val(p.y,j,0,i) - 0.5);
       alpha += n.y*a + n.z*(i - 0.5) + n.x*(j - 0.5);
       ni++;
     }
     
#line 320
if (val(p.z,i,j,0) > 0. && val(p.z,i,j,0) < 1.) {
       double a = sign(val(Phi,i,j,0) - val)*(val(p.z,i,j,0) - 0.5);
       alpha += n.z*a + n.x*(i - 0.5) + n.y*(j - 0.5);
       ni++;
     }}




      if (ni == 0)
 val(c,0,0,0) = val(as.x,0,0,0);
      else if (ni < 3 || ni > 6)
 val(c,0,0,0) = 0.;
      else
 val(c,0,0,0) = plane_volume (n, alpha/ni);
    }
  }end_foreach();}END_FOREACH delete((scalar*)((vector[]){p,{{-1},{-1},{-1}}}));
#line 351 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
end_tracing("fractions","/home/esteban/ProgramFile/basilisk/src/fractions.h",351);}
#line 395 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  coord n;
  double nn = 0.;
  if (!(3 == 2)) qassert ("/home/esteban/ProgramFile/basilisk/src/fractions.h", 399, "dimension == 2");
   {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  } 
#line 400
{
    n.y = (val(c,0,-1,1) + 2.*val(c,0,-1,0) + val(c,0,-1,-1) -
    val(c,0,+1,1) - 2.*val(c,0,+1,0) - val(c,0,+1,-1));
    nn += fabs(n.y);
  } 
#line 400
{
    n.z = (val(c,1,0,-1) + 2.*val(c,0,0,-1) + val(c,-1,0,-1) -
    val(c,1,0,+1) - 2.*val(c,0,0,+1) - val(c,-1,0,+1));
    nn += fabs(n.z);
  }

  if (nn > 0.)
    {
      n.x /= nn;
      
#line 408
n.y /= nn;
      
#line 408
n.z /= nn;}
  else
    n.x = 1.;
  return n;
}





coord facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if (s.x.i >= 0) {
    coord n;
    double nn = 0.;
     {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    } 
#line 423
{
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    } 
#line 423
{
      n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
      nn += fabs(n.z);
    }
    if (nn > 0.)
      {
 n.x /= nn;
 
#line 429
n.y /= nn;
 
#line 429
n.z /= nn;}
    else
      {
 n.x = 1./3;
 
#line 432
n.y = 1./3;
 
#line 432
n.z = 1./3;}
    return n;
  }
  return mycs (point, c);
}






#line 418
static void _stencil_facet_normal (Point point, scalar c, vector s)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if (s.x.i >= 0) {    
    
    
     { 
_stencil_val(s.x,0,0,0); _stencil_val(s.x,1,0,0);  
       
       
    
#line 426
} 
#line 423
{ 
_stencil_val(s.y,0,0,0); _stencil_val(s.y,0,1,0);  
       
       
    
#line 426
} 
#line 423
{ 
_stencil_val(s.z,0,0,0); _stencil_val(s.z,0,0,1);  
       
       
    
#line 426
}
      
   
       
    
       
  
    return ;
  } 
_stencil_mycs (point, c);
  
#line 435
return;
}
#line 445 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
     
void reconstruction (const scalar c, vector n, scalar alpha)
{tracing("reconstruction","/home/esteban/ProgramFile/basilisk/src/fractions.h",446);
  foreach_stencil(1,{0},) {





_stencil_val(c,0,0,0); _stencil_val(c,0,0,0);{ {
      _stencil_val_a(alpha,0,0,0);  
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 457
{_stencil_val_a(n.y,0,0,0);  }
 
#line 457
{_stencil_val_a(n.z,0,0,0);  }
    } 
{  






       _stencil_mycs (point, c);
      
 {_stencil_val_a(n.x,0,0,0);  }
 
#line 468
{_stencil_val_a(n.y,0,0,0);  }
 
#line 468
{_stencil_val_a(n.z,0,0,0);  }
_stencil_val(c,0,0,0);
      
#line 469
_stencil_val_a(alpha,0,0,0);    
    }}





          
    
  
#line 471
}end_foreach_stencil()
   BEGIN_FOREACH{
#line 448
foreach() {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      
 val(n.x,0,0,0) = 0.;
 
#line 457
val(n.y,0,0,0) = 0.;
 
#line 457
val(n.z,0,0,0) = 0.;
    }
    else {






      coord m = mycs (point, c);
      
 val(n.x,0,0,0) = m.x;
 
#line 468
val(n.y,0,0,0) = m.y;
 
#line 468
val(n.z,0,0,0) = m.z;
      val(alpha,0,0,0) = plane_alpha (val(c,0,0,0), m);
    }
  }end_foreach();}END_FOREACH 
#line 489 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
end_tracing("reconstruction","/home/esteban/ProgramFile/basilisk/src/fractions.h",489);}
#line 509 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
     
void output_facets (scalar c, FILE * fp, vector s)
{tracing("output_facets","/home/esteban/ProgramFile/basilisk/src/fractions.h",510);
  foreach_stencil (0,{0},)
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {  
       _stencil_facet_normal (point, c, s);     
      _stencil_val(c,0,0,0);        
#line 525 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
      
                 
              
     
 
        
 

    }        }end_foreach_stencil()
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
 BEGIN_FOREACH{
#line 512
foreach ()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = plane_alpha (val(c,0,0,0), n);
#line 525 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++)
 fprintf (fp, "%g %g %g\n",
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
 fputc ('\n', fp);

    }end_foreach();}END_FOREACH 
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif


  
#line 535
fflush (fp);
end_tracing("output_facets","/home/esteban/ProgramFile/basilisk/src/fractions.h",536);}







     
double interface_area (scalar c)
{tracing("interface_area","/home/esteban/ProgramFile/basilisk/src/fractions.h",545);
  double area = 0.;
  foreach_stencil (1,{0},)
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
       _stencil_mycs (point, c);     
      _stencil_val(c,0,0,0); 
             
    }        }end_foreach_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel reduction(+:area)){
#line 548
foreach ()
    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = plane_alpha (val(c,0,0,0), n);
      area += pow(Delta, 3 - 1)*plane_area_center (n, alpha, &p);
    }end_foreach();mpi_all_reduce_array(&area,double,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
  
#line 554
{end_tracing("interface_area","/home/esteban/ProgramFile/basilisk/src/fractions.h",554);return area;}
end_tracing("interface_area","/home/esteban/ProgramFile/basilisk/src/fractions.h",555);}
#line 36 "/home/esteban/ProgramFile/basilisk/src/vof.h"
#line 44 "/home/esteban/ProgramFile/basilisk/src/vof.h"
extern scalar * interfaces;
extern vector uf;
extern double dt;








static double vof_concentration_gradient_x (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  static const double cmin = 0.5;
  double cl = val(c,-1,0,0), cc = val(c,0,0,0), cr = val(c,1,0,0);
  if (_attribute[t.i].inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && _attribute[t.i].gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
 if (_attribute[t.i].gradient)
   return _attribute[t.i].gradient (val(t,-1,0,0)/cl, val(t,0,0,0)/cc, val(t,1,0,0)/cr)/Delta;
 else
   return (val(t,1,0,0)/cr - val(t,-1,0,0)/cl)/(2.*Delta);
      }
      else
 return (val(t,1,0,0)/cr - val(t,0,0,0)/cc)/Delta;
    }
    else if (cl >= cmin)
      return (val(t,0,0,0)/cc - val(t,-1,0,0)/cl)/Delta;
  }
  return 0.;
}

#line 55
static double vof_concentration_gradient_y (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  static const double cmin = 0.5;
  double cl = val(c,0,-1,0), cc = val(c,0,0,0), cr = val(c,0,1,0);
  if (_attribute[t.i].inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && _attribute[t.i].gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
 if (_attribute[t.i].gradient)
   return _attribute[t.i].gradient (val(t,0,-1,0)/cl, val(t,0,0,0)/cc, val(t,0,1,0)/cr)/Delta;
 else
   return (val(t,0,1,0)/cr - val(t,0,-1,0)/cl)/(2.*Delta);
      }
      else
 return (val(t,0,1,0)/cr - val(t,0,0,0)/cc)/Delta;
    }
    else if (cl >= cmin)
      return (val(t,0,0,0)/cc - val(t,0,-1,0)/cl)/Delta;
  }
  return 0.;
}

#line 55
static double vof_concentration_gradient_z (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  static const double cmin = 0.5;
  double cl = val(c,0,0,-1), cc = val(c,0,0,0), cr = val(c,0,0,1);
  if (_attribute[t.i].inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && _attribute[t.i].gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
 if (_attribute[t.i].gradient)
   return _attribute[t.i].gradient (val(t,0,0,-1)/cl, val(t,0,0,0)/cc, val(t,0,0,1)/cr)/Delta;
 else
   return (val(t,0,0,1)/cr - val(t,0,0,-1)/cl)/(2.*Delta);
      }
      else
 return (val(t,0,0,1)/cr - val(t,0,0,0)/cc)/Delta;
    }
    else if (cl >= cmin)
      return (val(t,0,0,0)/cc - val(t,0,0,-1)/cl)/Delta;
  }
  return 0.;
}









#line 55
static void _stencil_vof_concentration_gradient_x (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;           
  
   _stencil_val(c,1,0,0); _stencil_val(c,0,0,0); _stencil_val(c,-1,0,0); 


{
{ {
{ {
 if (_attribute[t.i].gradient)
   {_stencil_val(t,-1,0,0); _stencil_val(t,0,0,0); _stencil_val(t,1,0,0);  }
 else
   {_stencil_val(t,1,0,0); _stencil_val(t,-1,0,0);  }
      }
 
{_stencil_val(t,1,0,0); _stencil_val(t,0,0,0);  }}
         
      
    
#line 71
}
      
{_stencil_val(t,0,0,0); _stencil_val(t,-1,0,0);  }}
       
        
  
#line 74
} 
  
                  
         
  
#line 75
return ;
}

#line 55
static void _stencil_vof_concentration_gradient_y (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;           
  
   _stencil_val(c,0,1,0); _stencil_val(c,0,0,0); _stencil_val(c,0,-1,0); 


{
{ {
{ {
 if (_attribute[t.i].gradient)
   {_stencil_val(t,0,-1,0); _stencil_val(t,0,0,0); _stencil_val(t,0,1,0);  }
 else
   {_stencil_val(t,0,1,0); _stencil_val(t,0,-1,0);  }
      }
 
{_stencil_val(t,0,1,0); _stencil_val(t,0,0,0);  }}
         
      
    
#line 71
}
      
{_stencil_val(t,0,0,0); _stencil_val(t,0,-1,0);  }}
       
        
  
#line 74
} 
  
                  
         
  
#line 75
return ;
}

#line 55
static void _stencil_vof_concentration_gradient_z (Point point, scalar c, scalar t)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;           
  
   _stencil_val(c,0,0,1); _stencil_val(c,0,0,0); _stencil_val(c,0,0,-1); 


{
{ {
{ {
 if (_attribute[t.i].gradient)
   {_stencil_val(t,0,0,-1); _stencil_val(t,0,0,0); _stencil_val(t,0,0,1);  }
 else
   {_stencil_val(t,0,0,1); _stencil_val(t,0,0,-1);  }
      }
 
{_stencil_val(t,0,0,1); _stencil_val(t,0,0,0);  }}
         
      
    
#line 71
}
      
{_stencil_val(t,0,0,0); _stencil_val(t,0,0,-1);  }}
       
        
  
#line 74
} 
  
                  
         
  
#line 75
return ;
}
#line 127
static int defaults_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}
#line 127 "/home/esteban/ProgramFile/basilisk/src/vof.h"
      static int defaults_1(const int i,const double t,Event *_ev){tracing("defaults_1","/home/esteban/ProgramFile/basilisk/src/vof.h",127);
{
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar c=*_i;(&c)->i>=0;c=*++_i){ {
    scalar * tracers = _attribute[c.i].tracers;
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){
      _attribute[t.i].depends = list_add (_attribute[t.i].depends, c);}}
  }}}
}{end_tracing("defaults_1","/home/esteban/ProgramFile/basilisk/src/vof.h",134);return 0;}end_tracing("defaults_1","/home/esteban/ProgramFile/basilisk/src/vof.h",134);}





static int stability_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}






#line 140
      static int stability_0(const int i,const double t,Event *_ev){tracing("stability_0","/home/esteban/ProgramFile/basilisk/src/vof.h",140); {
  if (CFL > 0.5)
    CFL = 0.5;
}{end_tracing("stability_0","/home/esteban/ProgramFile/basilisk/src/vof.h",143);return 0;}end_tracing("stability_0","/home/esteban/ProgramFile/basilisk/src/vof.h",143);}
#line 157 "/home/esteban/ProgramFile/basilisk/src/vof.h"

static void sweep_x (scalar c, scalar cc, scalar * tcl)
{
  vector  n=new_vector("n");
  scalar  alpha=new_scalar("alpha"),  flux=new_scalar("flux");
  double cfl = 0.;
#line 171 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){ {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }}}




    foreach_stencil(1,{0},) {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 { _stencil_vof_concentration_gradient_x (point, c, t);_stencil_val_a(gf,0,0,0); }}}
    }end_foreach_stencil()




     BEGIN_FOREACH{
#line 182
foreach() {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 val(gf,0,0,0) = vof_concentration_gradient_x (point, c, t);}}
    }end_foreach();}END_FOREACH 
  }






  reconstruction (c, n, alpha);
  if(!is_constant(fm.x) && !is_constant(cm)){
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_x(){ {       






    _stencil_val(fm.x,0,0,0); _stencil_val(uf.x,0,0,0);     
    







_stencil_val(fm.x,0,0,0);_stencil_val(cm,0,0,0);
      {_stencil_val(fm.x,0,0,0);_stencil_val(cm,0,0,0);    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha, o_stencil,0,0);_stencil_val(n.z, o_stencil,0,0);_stencil_val(n.y, o_stencil,0,0);_stencil_val(n.x,o_stencil,0,0);
#line 225
_stencil_val(c, o_stencil,0,0);_stencil_val(c, o_stencil,0,0);_stencil_val(c,o_stencil,0,0);








_stencil_val(uf.x,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c, o_stencil,0,0);


{ {       
 _stencil_val(gf,o_stencil,0,0);_stencil_val(t, o_stencil,0,0);
_stencil_val(uf.x,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_x()end_foreach_face_stencil()
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*val(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val(fm.x,0,0,0)*s/(val(cm,0,0,0) + 0.) > cfl)
      cfl = un*val(fm.x,0,0,0)*s/(val(cm,0,0,0) + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}else if(is_constant(fm.x) && !is_constant(cm)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_x(){ {






; _stencil_val(uf.x,0,0,0);








;_stencil_val(cm,0,0,0);
      {;_stencil_val(cm,0,0,0);    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha, o_stencil,0,0);_stencil_val(n.z, o_stencil,0,0);_stencil_val(n.y, o_stencil,0,0);_stencil_val(n.x,o_stencil,0,0);
#line 225
_stencil_val(c, o_stencil,0,0);_stencil_val(c, o_stencil,0,0);_stencil_val(c,o_stencil,0,0);








_stencil_val(uf.x,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c, o_stencil,0,0);


{ {       
 _stencil_val(gf,o_stencil,0,0);_stencil_val(t, o_stencil,0,0);
_stencil_val(uf.x,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_x()end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*_const_fm.x + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*_const_fm.x*s/(val(cm,0,0,0) + 0.) > cfl)
      cfl = un*_const_fm.x*s/(val(cm,0,0,0) + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}else if(!is_constant(fm.x) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_x(){ {       






    _stencil_val(fm.x,0,0,0); _stencil_val(uf.x,0,0,0);     
    







_stencil_val(fm.x,0,0,0);;
      {_stencil_val(fm.x,0,0,0);;    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha, o_stencil,0,0);_stencil_val(n.z, o_stencil,0,0);_stencil_val(n.y, o_stencil,0,0);_stencil_val(n.x,o_stencil,0,0);
#line 225
_stencil_val(c, o_stencil,0,0);_stencil_val(c, o_stencil,0,0);_stencil_val(c,o_stencil,0,0);








_stencil_val(uf.x,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c, o_stencil,0,0);


{ {       
 _stencil_val(gf,o_stencil,0,0);_stencil_val(t, o_stencil,0,0);
_stencil_val(uf.x,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_x()end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*val(fm.x,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val(fm.x,0,0,0)*s/(_const_cm + 0.) > cfl)
      cfl = un*val(fm.x,0,0,0)*s/(_const_cm + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_x(){ {






; _stencil_val(uf.x,0,0,0);








;;
      {;;    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha, o_stencil,0,0);_stencil_val(n.z, o_stencil,0,0);_stencil_val(n.y, o_stencil,0,0);_stencil_val(n.x,o_stencil,0,0);
#line 225
_stencil_val(c, o_stencil,0,0);_stencil_val(c, o_stencil,0,0);_stencil_val(c,o_stencil,0,0);








_stencil_val(uf.x,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c, o_stencil,0,0);


{ {       
 _stencil_val(gf,o_stencil,0,0);_stencil_val(t, o_stencil,0,0);
_stencil_val(uf.x,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_x()end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_x(){ {






    double un = val(uf.x,0,0,0)*dt/(Delta*_const_fm.x + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*_const_fm.x*s/(_const_cm + 0.) > cfl)
      cfl = un*_const_fm.x*s/(_const_cm + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,i,0,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,i,0,0)/ci + s*min(1., 1. - s*un)*val(gf,i,0,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.x,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_x()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);




  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "src/vof.h:%d: warning: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      262, cfl - 0.5), fflush (ferr);
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  if(!is_constant(cm)){
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  foreach_stencil(1,{0},) {
_stencil_val(flux,0,0,0); _stencil_val(flux,1,0,0); _stencil_val(cc,0,0,0);_stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);_stencil_val(cm,0,0,0);
    
#line 284
_stencil_val_r(c,0,0,0);     





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      {_stencil_val(tflux,0,0,0); _stencil_val(tflux,1,0,0); _stencil_val(tc,0,0,0);_stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);_stencil_val(cm,0,0,0);_stencil_val_r(t,0,0,0);     }}}

  }end_foreach_stencil() BEGIN_FOREACH{
#line 283
foreach() {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val(cm,0,0,0)*Delta);





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0) + val(tc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val(cm,0,0,0)*Delta);}}

  }end_foreach();}END_FOREACH }else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  foreach_stencil(1,{0},) {
_stencil_val(flux,0,0,0); _stencil_val(flux,1,0,0); _stencil_val(cc,0,0,0);_stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);;
    
#line 284
_stencil_val_r(c,0,0,0);     





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      {_stencil_val(tflux,0,0,0); _stencil_val(tflux,1,0,0); _stencil_val(tc,0,0,0);_stencil_val(uf.x,1,0,0); _stencil_val(uf.x,0,0,0);;_stencil_val_r(t,0,0,0);     }}}

  }end_foreach_stencil()
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
   BEGIN_FOREACH{foreach() {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(_const_cm*Delta);





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,1,0,0) + val(tc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(_const_cm*Delta);}}

  }end_foreach();}END_FOREACH }
#line 317 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){flux,alpha,n.x,n.y,n.z,{-1}}));
}

#line 158
static void sweep_y (scalar c, scalar cc, scalar * tcl)
{
  vector  n=new_vector("n");
  scalar  alpha=new_scalar("alpha"),  flux=new_scalar("flux");
  double cfl = 0.;
#line 171 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){ {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }}}




    foreach_stencil(1,{0},) {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 { _stencil_vof_concentration_gradient_y (point, c, t);_stencil_val_a(gf,0,0,0); }}}
    }end_foreach_stencil()




     BEGIN_FOREACH{
#line 182
foreach() {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 val(gf,0,0,0) = vof_concentration_gradient_y (point, c, t);}}
    }end_foreach();}END_FOREACH 
  }






  reconstruction (c, n, alpha);
  if(!is_constant(fm.y) && !is_constant(cm)){
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_y(){ {       






    _stencil_val(fm.y,0,0,0); _stencil_val(uf.y,0,0,0);     
    







_stencil_val(fm.y,0,0,0);_stencil_val(cm,0,0,0);
      {_stencil_val(fm.y,0,0,0);_stencil_val(cm,0,0,0);    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha,0, o_stencil,0);_stencil_val(n.x,0, o_stencil,0);_stencil_val(n.z,0, o_stencil,0);_stencil_val(n.y,0,o_stencil,0);
#line 225
_stencil_val(c,0, o_stencil,0);_stencil_val(c,0, o_stencil,0);_stencil_val(c,0,o_stencil,0);








_stencil_val(uf.y,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0, o_stencil,0);


{ {       
 _stencil_val(gf,0,o_stencil,0);_stencil_val(t,0, o_stencil,0);
_stencil_val(uf.y,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_y()end_foreach_face_stencil()
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*val(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val(fm.y,0,0,0)*s/(val(cm,0,0,0) + 0.) > cfl)
      cfl = un*val(fm.y,0,0,0)*s/(val(cm,0,0,0) + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}else if(is_constant(fm.y) && !is_constant(cm)){_coord _const_fm={_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX],_constant[fm.x.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_y(){ {






; _stencil_val(uf.y,0,0,0);








;_stencil_val(cm,0,0,0);
      {;_stencil_val(cm,0,0,0);    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha,0, o_stencil,0);_stencil_val(n.x,0, o_stencil,0);_stencil_val(n.z,0, o_stencil,0);_stencil_val(n.y,0,o_stencil,0);
#line 225
_stencil_val(c,0, o_stencil,0);_stencil_val(c,0, o_stencil,0);_stencil_val(c,0,o_stencil,0);








_stencil_val(uf.y,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0, o_stencil,0);


{ {       
 _stencil_val(gf,0,o_stencil,0);_stencil_val(t,0, o_stencil,0);
_stencil_val(uf.y,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_y()end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*_const_fm.y + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*_const_fm.y*s/(val(cm,0,0,0) + 0.) > cfl)
      cfl = un*_const_fm.y*s/(val(cm,0,0,0) + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}else if(!is_constant(fm.y) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_y(){ {       






    _stencil_val(fm.y,0,0,0); _stencil_val(uf.y,0,0,0);     
    







_stencil_val(fm.y,0,0,0);;
      {_stencil_val(fm.y,0,0,0);;    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha,0, o_stencil,0);_stencil_val(n.x,0, o_stencil,0);_stencil_val(n.z,0, o_stencil,0);_stencil_val(n.y,0,o_stencil,0);
#line 225
_stencil_val(c,0, o_stencil,0);_stencil_val(c,0, o_stencil,0);_stencil_val(c,0,o_stencil,0);








_stencil_val(uf.y,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0, o_stencil,0);


{ {       
 _stencil_val(gf,0,o_stencil,0);_stencil_val(t,0, o_stencil,0);
_stencil_val(uf.y,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_y()end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*val(fm.y,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val(fm.y,0,0,0)*s/(_const_cm + 0.) > cfl)
      cfl = un*val(fm.y,0,0,0)*s/(_const_cm + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}else {_coord _const_fm={_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX],_constant[fm.x.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_y(){ {






; _stencil_val(uf.y,0,0,0);








;;
      {;;    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha,0, o_stencil,0);_stencil_val(n.x,0, o_stencil,0);_stencil_val(n.z,0, o_stencil,0);_stencil_val(n.y,0,o_stencil,0);
#line 225
_stencil_val(c,0, o_stencil,0);_stencil_val(c,0, o_stencil,0);_stencil_val(c,0,o_stencil,0);








_stencil_val(uf.y,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0, o_stencil,0);


{ {       
 _stencil_val(gf,0,o_stencil,0);_stencil_val(t,0, o_stencil,0);
_stencil_val(uf.y,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_y()end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_y(){ {






    double un = val(uf.y,0,0,0)*dt/(Delta*_const_fm.y + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*_const_fm.y*s/(_const_cm + 0.) > cfl)
      cfl = un*_const_fm.y*s/(_const_cm + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,i,0);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,i,0)/ci + s*min(1., 1. - s*un)*val(gf,0,i,0)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.y,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_y()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);




  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "src/vof.h:%d: warning: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      262, cfl - 0.5), fflush (ferr);
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  if(!is_constant(cm)){
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  foreach_stencil(1,{0},) {
_stencil_val(flux,0,0,0); _stencil_val(flux,0,1,0); _stencil_val(cc,0,0,0);_stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);_stencil_val(cm,0,0,0);
    
#line 284
_stencil_val_r(c,0,0,0);     





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      {_stencil_val(tflux,0,0,0); _stencil_val(tflux,0,1,0); _stencil_val(tc,0,0,0);_stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);_stencil_val(cm,0,0,0);_stencil_val_r(t,0,0,0);     }}}

  }end_foreach_stencil() BEGIN_FOREACH{
#line 283
foreach() {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val(cm,0,0,0)*Delta);





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0) + val(tc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val(cm,0,0,0)*Delta);}}

  }end_foreach();}END_FOREACH }else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  foreach_stencil(1,{0},) {
_stencil_val(flux,0,0,0); _stencil_val(flux,0,1,0); _stencil_val(cc,0,0,0);_stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);;
    
#line 284
_stencil_val_r(c,0,0,0);     





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      {_stencil_val(tflux,0,0,0); _stencil_val(tflux,0,1,0); _stencil_val(tc,0,0,0);_stencil_val(uf.y,0,1,0); _stencil_val(uf.y,0,0,0);;_stencil_val_r(t,0,0,0);     }}}

  }end_foreach_stencil()
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
   BEGIN_FOREACH{foreach() {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(_const_cm*Delta);





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,1,0) + val(tc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(_const_cm*Delta);}}

  }end_foreach();}END_FOREACH }
#line 317 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){flux,alpha,n.x,n.y,n.z,{-1}}));
}

#line 158
static void sweep_z (scalar c, scalar cc, scalar * tcl)
{
  vector  n=new_vector("n");
  scalar  alpha=new_scalar("alpha"),  flux=new_scalar("flux");
  double cfl = 0.;
#line 171 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  scalar * tracers = _attribute[c.i].tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){ {
      scalar gf = new_scalar("gf"), flux = new_scalar("flux");
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }}}




    foreach_stencil(1,{0},) {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 { _stencil_vof_concentration_gradient_z (point, c, t);_stencil_val_a(gf,0,0,0); }}}
    }end_foreach_stencil()




     BEGIN_FOREACH{
#line 182
foreach() {
      scalar t, gf;
      {scalar*_i0=gfl;scalar*_i1= tracers;if(_i0)for(gf=*_i0,t=*_i1;_i0->i>= 0;gf=*++_i0,t=*++_i1){
 val(gf,0,0,0) = vof_concentration_gradient_z (point, c, t);}}
    }end_foreach();}END_FOREACH 
  }






  reconstruction (c, n, alpha);
  if(!is_constant(fm.z) && !is_constant(cm)){
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_z(){ {       






    _stencil_val(fm.z,0,0,0); _stencil_val(uf.z,0,0,0);     
    







_stencil_val(fm.z,0,0,0);_stencil_val(cm,0,0,0);
      {_stencil_val(fm.z,0,0,0);_stencil_val(cm,0,0,0);    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha,0,0, o_stencil);_stencil_val(n.y,0,0, o_stencil);_stencil_val(n.x,0,0, o_stencil);_stencil_val(n.z,0,0,o_stencil);
#line 225
_stencil_val(c,0,0, o_stencil);_stencil_val(c,0,0, o_stencil);_stencil_val(c,0,0,o_stencil);








_stencil_val(uf.z,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0,0, o_stencil);


{ {       
 _stencil_val(gf,0,0,o_stencil);_stencil_val(t,0,0, o_stencil);
_stencil_val(uf.z,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_z()end_foreach_face_stencil()
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_z(){ {






    double un = val(uf.z,0,0,0)*dt/(Delta*val(fm.z,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val(fm.z,0,0,0)*s/(val(cm,0,0,0) + 0.) > cfl)
      cfl = un*val(fm.z,0,0,0)*s/(val(cm,0,0,0) + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,0,i);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,0,i)/ci + s*min(1., 1. - s*un)*val(gf,0,0,i)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.z,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_z()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}else if(is_constant(fm.z) && !is_constant(cm)){_coord _const_fm={_constant[fm.z.i-_NVARMAX],_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_z(){ {






; _stencil_val(uf.z,0,0,0);








;_stencil_val(cm,0,0,0);
      {;_stencil_val(cm,0,0,0);    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha,0,0, o_stencil);_stencil_val(n.y,0,0, o_stencil);_stencil_val(n.x,0,0, o_stencil);_stencil_val(n.z,0,0,o_stencil);
#line 225
_stencil_val(c,0,0, o_stencil);_stencil_val(c,0,0, o_stencil);_stencil_val(c,0,0,o_stencil);








_stencil_val(uf.z,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0,0, o_stencil);


{ {       
 _stencil_val(gf,0,0,o_stencil);_stencil_val(t,0,0, o_stencil);
_stencil_val(uf.z,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_z()end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_z(){ {






    double un = val(uf.z,0,0,0)*dt/(Delta*_const_fm.z + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*_const_fm.z*s/(val(cm,0,0,0) + 0.) > cfl)
      cfl = un*_const_fm.z*s/(val(cm,0,0,0) + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,0,i);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,0,i)/ci + s*min(1., 1. - s*un)*val(gf,0,0,i)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.z,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_z()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}else if(!is_constant(fm.z) && is_constant(cm)){double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_z(){ {       






    _stencil_val(fm.z,0,0,0); _stencil_val(uf.z,0,0,0);     
    







_stencil_val(fm.z,0,0,0);;
      {_stencil_val(fm.z,0,0,0);;    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha,0,0, o_stencil);_stencil_val(n.y,0,0, o_stencil);_stencil_val(n.x,0,0, o_stencil);_stencil_val(n.z,0,0,o_stencil);
#line 225
_stencil_val(c,0,0, o_stencil);_stencil_val(c,0,0, o_stencil);_stencil_val(c,0,0,o_stencil);








_stencil_val(uf.z,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0,0, o_stencil);


{ {       
 _stencil_val(gf,0,0,o_stencil);_stencil_val(t,0,0, o_stencil);
_stencil_val(uf.z,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_z()end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_z(){ {






    double un = val(uf.z,0,0,0)*dt/(Delta*val(fm.z,0,0,0) + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*val(fm.z,0,0,0)*s/(_const_cm + 0.) > cfl)
      cfl = un*val(fm.z,0,0,0)*s/(_const_cm + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,0,i);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,0,i)/ci + s*min(1., 1. - s*un)*val(gf,0,0,i)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.z,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_z()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}else {_coord _const_fm={_constant[fm.z.i-_NVARMAX],_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX]};NOT_UNUSED(_const_fm);double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
  
#line 195
foreach_face_stencil(1,{0},)_stencil_is_face_z(){ {






; _stencil_val(uf.z,0,0,0);








;;
      {;;    }              
       
      
      







         
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    
_stencil_val(alpha,0,0, o_stencil);_stencil_val(n.y,0,0, o_stencil);_stencil_val(n.x,0,0, o_stencil);_stencil_val(n.z,0,0,o_stencil);
#line 225
_stencil_val(c,0,0, o_stencil);_stencil_val(c,0,0, o_stencil);_stencil_val(c,0,0,o_stencil);








_stencil_val(uf.z,0,0,0);





    
#line 234
_stencil_val_a(flux,0,0,0);  






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {     
      _stencil_val(c,0,0, o_stencil);


{ {       
 _stencil_val(gf,0,0,o_stencil);_stencil_val(t,0,0, o_stencil);
_stencil_val(uf.z,0,0,0);
 
#line 248
_stencil_val_a(tflux,0,0,0);  
      }
 
{_stencil_val_a(tflux,0,0,0);  }} 
      
          
         
      
    
#line 252
}}}
  }}end__stencil_is_face_z()end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction (max:cfl)){
#line 195
foreach_face_generic()is_face_z(){ {






    double un = val(uf.z,0,0,0)*dt/(Delta*_const_fm.z + 0.), s = sign(un);
    int i = -(s + 1.)/2.;







    if (un*_const_fm.z*s/(_const_cm + 0.) > cfl)
      cfl = un*_const_fm.z*s/(_const_cm + 0.);
#line 225 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);






    scalar t, gf, tflux;
    {scalar*_i0=tfluxl;scalar*_i1=gfl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,gf=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,gf=*++_i1,t=*++_i2){ {
      double cf1 = cf, ci = val(c,0,0,i);
      if (_attribute[t.i].inverse)
 cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
 double ff = val(t,0,0,i)/ci + s*min(1., 1. - s*un)*val(gf,0,0,i)*Delta/2.;
 val(tflux,0,0,0) = ff*cf1*val(uf.z,0,0,0);
      }
      else
 val(tflux,0,0,0) = 0.;
    }}}
  }}end_is_face_z()end_foreach_face_generic();mpi_all_reduce_array(&cfl,double,MPI_MAX,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 253
}
  delete (gfl); pfree (gfl,__func__,__FILE__,__LINE__);




  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "src/vof.h:%d: warning: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      262, cfl - 0.5), fflush (ferr);
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  if(!is_constant(cm)){
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  foreach_stencil(1,{0},) {
_stencil_val(flux,0,0,0); _stencil_val(flux,0,0,1); _stencil_val(cc,0,0,0);_stencil_val(uf.z,0,0,1); _stencil_val(uf.z,0,0,0);_stencil_val(cm,0,0,0);
    
#line 284
_stencil_val_r(c,0,0,0);     





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      {_stencil_val(tflux,0,0,0); _stencil_val(tflux,0,0,1); _stencil_val(tc,0,0,0);_stencil_val(uf.z,0,0,1); _stencil_val(uf.z,0,0,0);_stencil_val(cm,0,0,0);_stencil_val_r(t,0,0,0);     }}}

  }end_foreach_stencil() BEGIN_FOREACH{
#line 283
foreach() {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,0,1) + val(cc,0,0,0)*(val(uf.z,0,0,1) - val(uf.z,0,0,0)))/(val(cm,0,0,0)*Delta);





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,0,1) + val(tc,0,0,0)*(val(uf.z,0,0,1) - val(uf.z,0,0,0)))/(val(cm,0,0,0)*Delta);}}

  }end_foreach();}END_FOREACH }else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  foreach_stencil(1,{0},) {
_stencil_val(flux,0,0,0); _stencil_val(flux,0,0,1); _stencil_val(cc,0,0,0);_stencil_val(uf.z,0,0,1); _stencil_val(uf.z,0,0,0);;
    
#line 284
_stencil_val_r(c,0,0,0);     





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      {_stencil_val(tflux,0,0,0); _stencil_val(tflux,0,0,1); _stencil_val(tc,0,0,0);_stencil_val(uf.z,0,0,1); _stencil_val(uf.z,0,0,0);;_stencil_val_r(t,0,0,0);     }}}

  }end_foreach_stencil()
#line 283 "/home/esteban/ProgramFile/basilisk/src/vof.h"
   BEGIN_FOREACH{foreach() {
    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,0,1) + val(cc,0,0,0)*(val(uf.z,0,0,1) - val(uf.z,0,0,0)))/(_const_cm*Delta);





    scalar t, tc, tflux;
    {scalar*_i0= tfluxl;scalar*_i1= tcl;scalar*_i2= tracers;if(_i0)for(tflux=*_i0,tc=*_i1,t=*_i2;_i0->i>= 0;tflux=*++_i0,tc=*++_i1,t=*++_i2){
      val(t,0,0,0) += dt*(val(tflux,0,0,0) - val(tflux,0,0,1) + val(tc,0,0,0)*(val(uf.z,0,0,1) - val(uf.z,0,0,0)))/(_const_cm*Delta);}}

  }end_foreach();}END_FOREACH }
#line 317 "/home/esteban/ProgramFile/basilisk/src/vof.h"
  delete (tfluxl); pfree (tfluxl,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){flux,alpha,n.x,n.y,n.z,{-1}}));
}






void vof_advection (scalar * interfaces, int i)
{
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar c=*_i;(&c)->i>=0;c=*++_i){ {
#line 337 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    scalar  cc=new_scalar("cc"), * tcl = NULL, * tracers = _attribute[c.i].tracers;
    {scalar*_i=(scalar*)( tracers);if(_i)for(scalar t=*_i;(&t)->i>=0;t=*++_i){ {

      scalar tc = new_scalar("tc");
      tcl = list_append (tcl, tc);
#line 351 "/home/esteban/ProgramFile/basilisk/src/vof.h"
    }}}
    foreach_stencil(1,{0},) {
_stencil_val(c,0,0,0);
      
#line 353
_stencil_val_a(cc,0,0,0);    

      scalar t, tc;
      {scalar*_i0= tcl;scalar*_i1= tracers;if(_i0)for(tc=*_i0,t=*_i1;_i0->i>= 0;tc=*++_i0,t=*++_i1){ {
 if (_attribute[t.i].inverse)
   { _stencil_val(c,0,0,0); _stencil_val(t,0,0,0); _stencil_val(c,0,0,0);_stencil_val_a(tc,0,0,0);       }
 else
   { _stencil_val(c,0,0,0); _stencil_val(t,0,0,0);_stencil_val(c,0,0,0);_stencil_val_a(tc,0,0,0);      }
      }}}

    }end_foreach_stencil()
     BEGIN_FOREACH{
#line 352
foreach() {
      val(cc,0,0,0) = (val(c,0,0,0) > 0.5);

      scalar t, tc;
      {scalar*_i0= tcl;scalar*_i1= tracers;if(_i0)for(tc=*_i0,t=*_i1;_i0->i>= 0;tc=*++_i0,t=*++_i1){ {
 if (_attribute[t.i].inverse)
   val(tc,0,0,0) = val(c,0,0,0) < 0.5 ? val(t,0,0,0)/(1. - val(c,0,0,0)) : 0.;
 else
   val(tc,0,0,0) = val(c,0,0,0) > 0.5 ? val(t,0,0,0)/val(c,0,0,0) : 0.;
      }}}

    }end_foreach();}END_FOREACH 






    void (* sweep[3]) (scalar, scalar, scalar *);
    int d = 0;
    
      sweep[d++] = sweep_x;
      
#line 373
sweep[d++] = sweep_y;
      
#line 373
sweep[d++] = sweep_z;
    for (d = 0; d < 3; d++)
      sweep[(i + d) % 3] (c, cc, tcl);
    delete (tcl), pfree (tcl,__func__,__FILE__,__LINE__);delete((scalar*)((scalar[]){cc,{-1}}));
  }}}
}

static int vof_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 380
      static int vof_0(const int i,const double t,Event *_ev){tracing("vof_0","/home/esteban/ProgramFile/basilisk/src/vof.h",380);
{vof_advection (interfaces, i);
  
#line 381
}{end_tracing("vof_0","/home/esteban/ProgramFile/basilisk/src/vof.h",381);return 0;}end_tracing("vof_0","/home/esteban/ProgramFile/basilisk/src/vof.h",381);}
#line 14 "/home/esteban/ProgramFile/basilisk/src/two-phase.h"

scalar  f={11}, * interfaces =((scalar[]) {{11},{-1}});

#line 1 "two-phase-generic.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h"
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;





vector  alphav={{12},{13},{14}};
scalar  rhov={15};

static int defaults_2_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}


#line 10
      static int defaults_2(const int i,const double t,Event *_ev){tracing("defaults_2","/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h",10);
{
  alpha = alphav;
  rho = rhov;





  if (mu1 || mu2)
    mu = new_face_vector("mu");




  display ("draw_vof (c = 'f');"
#line 494 "/home/esteban/ProgramFile/basilisk/src/common.h"
, false
#line 25 "/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h"
);
}{end_tracing("defaults_2","/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h",26);return 0;}end_tracing("defaults_2","/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h",26);}
#line 50
static int tracer_advection_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 50 "/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h"
static int tracer_advection_0(const int i,const double t,Event *_ev){
{
#line 79 "/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h"
}return 0;}



static int properties_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}




#line 83
      static int properties_0(const int i,const double t,Event *_ev){tracing("properties_0","/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h",83);
{
  if(!is_constant(fm.x)){
  
#line 85
foreach_face_stencil(1,{0},){_stencil_is_face_x(){ {    
     _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0); 
_stencil_val(fm.x,0,0,0);
    
#line 87
_stencil_val_a(alphav.x,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu; 
_stencil_val(fm.x,0,0,0);
      
#line 90
_stencil_val_a(muv.x,0,0,0);     
    }
  }}end__stencil_is_face_x()
#line 85
_stencil_is_face_y(){ {    
     _stencil_val(f,0,-1,0);_stencil_val(f,0,0,0); 
_stencil_val(fm.y,0,0,0);
    
#line 87
_stencil_val_a(alphav.y,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu; 
_stencil_val(fm.y,0,0,0);
      
#line 90
_stencil_val_a(muv.y,0,0,0);     
    }
  }}end__stencil_is_face_y()
#line 85
_stencil_is_face_z(){ {    
     _stencil_val(f,0,0,-1);_stencil_val(f,0,0,0); 
_stencil_val(fm.z,0,0,0);
    
#line 87
_stencil_val_a(alphav.z,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu; 
_stencil_val(fm.z,0,0,0);
      
#line 90
_stencil_val_a(muv.z,0,0,0);     
    }
  }}end__stencil_is_face_z()}end_foreach_face_stencil() BEGIN_FOREACH{
#line 85
foreach_face_generic(){is_face_x(){ {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = val(fm.x,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = val(fm.x,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_x()
#line 85
is_face_y(){ {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = val(fm.y,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = val(fm.y,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_y()
#line 85
is_face_z(){ {
    double ff = (val(f,0,0,0) + val(f,0,0,-1))/2.;
    val(alphav.z,0,0,0) = val(fm.z,0,0,0)/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.z,0,0,0) = val(fm.z,0,0,0)*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 85
foreach_face_stencil(1,{0},){_stencil_is_face_x(){ {    
     _stencil_val(f,-1,0,0);_stencil_val(f,0,0,0);
;
    
#line 87
_stencil_val_a(alphav.x,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu;
;
      
#line 90
_stencil_val_a(muv.x,0,0,0);     
    }
  }}end__stencil_is_face_x()
#line 85
_stencil_is_face_y(){ {    
     _stencil_val(f,0,-1,0);_stencil_val(f,0,0,0);
;
    
#line 87
_stencil_val_a(alphav.y,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu;
;
      
#line 90
_stencil_val_a(muv.y,0,0,0);     
    }
  }}end__stencil_is_face_y()
#line 85
_stencil_is_face_z(){ {    
     _stencil_val(f,0,0,-1);_stencil_val(f,0,0,0);
;
    
#line 87
_stencil_val_a(alphav.z,0,0,0);     
    if (mu1 || mu2) {
      vector muv = mu;
;
      
#line 90
_stencil_val_a(muv.z,0,0,0);     
    }
  }}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 85
foreach_face_generic(){is_face_x(){ {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = _const_fm.x/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.x,0,0,0) = _const_fm.x*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_x()
#line 85
is_face_y(){ {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = _const_fm.y/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.y,0,0,0) = _const_fm.y*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_y()
#line 85
is_face_z(){ {
    double ff = (val(f,0,0,0) + val(f,0,0,-1))/2.;
    val(alphav.z,0,0,0) = _const_fm.z/(clamp(ff,0.,1.)*(rho1 - rho2) + rho2);
    if (mu1 || mu2) {
      vector muv = mu;
      val(muv.z,0,0,0) = _const_fm.z*(clamp(ff,0.,1.)*(mu1 - mu2) + mu2);
    }
  }}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }

  if(!is_constant(cm)){

  
#line 94
foreach_stencil(1,{0},)
    { _stencil_val(cm,0,0,0);_stencil_val(f,0,0,0);_stencil_val_a(rhov,0,0,0);     }end_foreach_stencil() BEGIN_FOREACH{
#line 94
foreach()
    val(rhov,0,0,0) = val(cm,0,0,0)*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);end_foreach();}END_FOREACH }else {double _const_cm=_constant[cm.i-_NVARMAX];NOT_UNUSED(_const_cm);

  
#line 94
foreach_stencil(1,{0},)
    {;_stencil_val(f,0,0,0);_stencil_val_a(rhov,0,0,0);     }end_foreach_stencil()

   BEGIN_FOREACH{
#line 94
foreach()
    val(rhov,0,0,0) = _const_cm*(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);end_foreach();}END_FOREACH }





}{end_tracing("properties_0","/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h",101);return 0;}end_tracing("properties_0","/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h",101);}
#line 30 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
static int defaults_3_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}
#line 18 "/home/esteban/ProgramFile/basilisk/src/two-phase.h"
#line 14 "main3D.c"
#line 1 "tension.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/tension.h"
#line 15 "/home/esteban/ProgramFile/basilisk/src/tension.h"
#line 1 "iforce.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
#line 20 "/home/esteban/ProgramFile/basilisk/src/iforce.h"










      static int defaults_3(const int i,const double t,Event *_ev){tracing("defaults_3","/home/esteban/ProgramFile/basilisk/src/iforce.h",30); {
  if (is_constant(a.x)) {
    a = new_face_vector("a");
    if(!is_constant(a.x)){
    
#line 33
foreach_face_stencil(1,{0},){_stencil_is_face_x(){ {
      _stencil_val_a(a.x,0,0,0);
_stencil_val(a.x,0,0,0);     
      
    
#line 36
}}end__stencil_is_face_x()
#line 33
_stencil_is_face_y(){ {
      _stencil_val_a(a.y,0,0,0);
_stencil_val(a.y,0,0,0);     
      
    
#line 36
}}end__stencil_is_face_y()
#line 33
_stencil_is_face_z(){ {
      _stencil_val_a(a.z,0,0,0);
_stencil_val(a.z,0,0,0);     
      
    
#line 36
}}end__stencil_is_face_z()}end_foreach_face_stencil() BEGIN_FOREACH{
#line 33
foreach_face_generic(){is_face_x(){ {
      val(a.x,0,0,0) = 0.;
      dimensional (val(a.x,0,0,0) == Delta/sq(DT));
    }}end_is_face_x()
#line 33
is_face_y(){ {
      val(a.y,0,0,0) = 0.;
      dimensional (val(a.y,0,0,0) == Delta/sq(DT));
    }}end_is_face_y()
#line 33
is_face_z(){ {
      val(a.z,0,0,0) = 0.;
      dimensional (val(a.z,0,0,0) == Delta/sq(DT));
    }}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else {_coord _const_a={_constant[a.x.i-_NVARMAX],_constant[a.y.i-_NVARMAX],_constant[a.z.i-_NVARMAX]};NOT_UNUSED(_const_a);
    
#line 33
foreach_face_stencil(1,{0},){_stencil_is_face_x(){ {
      _stencil_val_a(a.x,0,0,0);
;     
      
    
#line 36
}}end__stencil_is_face_x()
#line 33
_stencil_is_face_y(){ {
      _stencil_val_a(a.y,0,0,0);
;     
      
    
#line 36
}}end__stencil_is_face_y()
#line 33
_stencil_is_face_z(){ {
      _stencil_val_a(a.z,0,0,0);
;     
      
    
#line 36
}}end__stencil_is_face_z()}end_foreach_face_stencil()
     BEGIN_FOREACH{
#line 33
foreach_face_generic(){is_face_x(){ {
      _const_a.x = 0.;
      dimensional (_const_a.x == Delta/sq(DT));
    }}end_is_face_x()
#line 33
is_face_y(){ {
      _const_a.y = 0.;
      dimensional (_const_a.y == Delta/sq(DT));
    }}end_is_face_y()
#line 33
is_face_z(){ {
      _const_a.z = 0.;
      dimensional (_const_a.z == Delta/sq(DT));
    }}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }
  }
}{end_tracing("defaults_3","/home/esteban/ProgramFile/basilisk/src/iforce.h",38);return 0;}end_tracing("defaults_3","/home/esteban/ProgramFile/basilisk/src/iforce.h",38);}






static int acceleration_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}







#line 45
      static int acceleration_0(const int i,const double t,Event *_ev){tracing("acceleration_0","/home/esteban/ProgramFile/basilisk/src/iforce.h",45);
{





  scalar * list = NULL;
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
    if (_attribute[f.i].phi.i) {
      list = list_add (list, f);






      foreach_stencil(1,{0},)
 {_stencil_val(f,0,0,0);_stencil_val_a(f,0,0,0);     }end_foreach_stencil()






       BEGIN_FOREACH{
#line 62
foreach()
 val(f,0,0,0) = clamp (val(f,0,0,0), 0., 1.);end_foreach();}END_FOREACH 
    }}}
#line 88 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
  vector ia = a;
  if(!is_constant(fm.x) && !is_constant(alpha.x)){
  
#line 89
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0); _stencil_val(fm.x,0,0,0); {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,-1,0,0);
   
#line 106
_stencil_val(phi,-1,0,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 





_stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);

 
#line 109
_stencil_val_r(ia.x,0,0,0);    
      }     }}}}end__stencil_is_face_x()
#line 89
_stencil_is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0); _stencil_val(fm.y,0,0,0); {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,-1,0);
   
#line 106
_stencil_val(phi,0,-1,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 





_stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);

 
#line 109
_stencil_val_r(ia.y,0,0,0);    
      }     }}}}end__stencil_is_face_y()
#line 89
_stencil_is_face_z(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,0,-1); _stencil_val(fm.z,0,0,0); {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,0,-1);
   
#line 106
_stencil_val(phi,0,0,-1); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,0,-1);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,0,-1);_stencil_val(phi,0,0,0); 





_stencil_val(alpha.z,0,0,0);_stencil_val(fm.z,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,0,0,-1);

 
#line 109
_stencil_val_r(ia.z,0,0,0);    
      }     }}}}end__stencil_is_face_z()}end_foreach_face_stencil() BEGIN_FOREACH{
#line 89
foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && val(fm.x,0,0,0) > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val(alpha.x,0,0,0)/(val(fm.x,0,0,0) + 0.)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 89
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && val(fm.y,0,0,0) > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val(alpha.y,0,0,0)/(val(fm.y,0,0,0) + 0.)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()
#line 89
is_face_z(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,0,-1) && val(fm.z,0,0,0) > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,0,-1) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,0,-1))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,0,-1) < 1e30 ? val(phi,0,0,-1) :
   0.;

 val(ia.z,0,0,0) += val(alpha.z,0,0,0)/(val(fm.z,0,0,0) + 0.)*phif*(val(f,0,0,0) - val(f,0,0,-1))/Delta;
      }}}}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 89
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);; {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,-1,0,0);
   
#line 106
_stencil_val(phi,-1,0,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 





_stencil_val(alpha.x,0,0,0);;_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);

 
#line 109
_stencil_val_r(ia.x,0,0,0);    
      }     }}}}end__stencil_is_face_x()
#line 89
_stencil_is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);; {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,-1,0);
   
#line 106
_stencil_val(phi,0,-1,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 





_stencil_val(alpha.y,0,0,0);;_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);

 
#line 109
_stencil_val_r(ia.y,0,0,0);    
      }     }}}}end__stencil_is_face_y()
#line 89
_stencil_is_face_z(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,0,-1);; {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,0,-1);
   
#line 106
_stencil_val(phi,0,0,-1); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,0,-1);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,0,-1);_stencil_val(phi,0,0,0); 





_stencil_val(alpha.z,0,0,0);;_stencil_val(f,0,0,0); _stencil_val(f,0,0,-1);

 
#line 109
_stencil_val_r(ia.z,0,0,0);    
      }     }}}}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 89
foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && _const_fm.x > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += val(alpha.x,0,0,0)/(_const_fm.x + 0.)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 89
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && _const_fm.y > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += val(alpha.y,0,0,0)/(_const_fm.y + 0.)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()
#line 89
is_face_z(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,0,-1) && _const_fm.z > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,0,-1) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,0,-1))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,0,-1) < 1e30 ? val(phi,0,0,-1) :
   0.;

 val(ia.z,0,0,0) += val(alpha.z,0,0,0)/(_const_fm.z + 0.)*phif*(val(f,0,0,0) - val(f,0,0,-1))/Delta;
      }}}}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else if(!is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 89
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0); _stencil_val(fm.x,0,0,0); {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,-1,0,0);
   
#line 106
_stencil_val(phi,-1,0,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0);





;_stencil_val(fm.x,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);

 
#line 109
_stencil_val_r(ia.x,0,0,0);    
      }     }}}}end__stencil_is_face_x()
#line 89
_stencil_is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0); _stencil_val(fm.y,0,0,0); {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,-1,0);
   
#line 106
_stencil_val(phi,0,-1,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0);





;_stencil_val(fm.y,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);

 
#line 109
_stencil_val_r(ia.y,0,0,0);    
      }     }}}}end__stencil_is_face_y()
#line 89
_stencil_is_face_z(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,0,-1); _stencil_val(fm.z,0,0,0); {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,0,-1);
   
#line 106
_stencil_val(phi,0,0,-1); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,0,-1);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,0,-1);_stencil_val(phi,0,0,0);





;_stencil_val(fm.z,0,0,0);_stencil_val(f,0,0,0); _stencil_val(f,0,0,-1);

 
#line 109
_stencil_val_r(ia.z,0,0,0);    
      }     }}}}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 89
foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && val(fm.x,0,0,0) > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += _const_alpha.x/(val(fm.x,0,0,0) + 0.)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 89
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && val(fm.y,0,0,0) > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += _const_alpha.y/(val(fm.y,0,0,0) + 0.)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()
#line 89
is_face_z(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,0,-1) && val(fm.z,0,0,0) > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,0,-1) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,0,-1))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,0,-1) < 1e30 ? val(phi,0,0,-1) :
   0.;

 val(ia.z,0,0,0) += _const_alpha.z/(val(fm.z,0,0,0) + 0.)*phif*(val(f,0,0,0) - val(f,0,0,-1))/Delta;
      }}}}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 89
foreach_face_stencil(1,{0},){_stencil_is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);; {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,-1,0,0);
   
#line 106
_stencil_val(phi,-1,0,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,-1,0,0);_stencil_val(phi,0,0,0);





;;_stencil_val(f,0,0,0); _stencil_val(f,-1,0,0);

 
#line 109
_stencil_val_r(ia.x,0,0,0);    
      }     }}}}end__stencil_is_face_x()
#line 89
_stencil_is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);; {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,-1,0);
   
#line 106
_stencil_val(phi,0,-1,0); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,-1,0);_stencil_val(phi,0,0,0);





;;_stencil_val(f,0,0,0); _stencil_val(f,0,-1,0);

 
#line 109
_stencil_val_r(ia.y,0,0,0);    
      }     }}}}end__stencil_is_face_y()
#line 89
_stencil_is_face_z(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      {_stencil_val(f,0,0,0); _stencil_val(f,0,0,-1);; {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;   
         
        
       
  
_stencil_val(phi,0,0,-1);
   
#line 106
_stencil_val(phi,0,0,-1); 
#line 105
_stencil_val(phi,0,0,0);
   
#line 105
_stencil_val(phi,0,0,0); 
#line 104
_stencil_val(phi,0,0,-1);_stencil_val(phi,0,0,0); 
#line 103
_stencil_val(phi,0,0,-1);_stencil_val(phi,0,0,0);





;;_stencil_val(f,0,0,0); _stencil_val(f,0,0,-1);

 
#line 109
_stencil_val_r(ia.z,0,0,0);    
      }     }}}}end__stencil_is_face_z()}end_foreach_face_stencil()
   BEGIN_FOREACH{
#line 89
foreach_face_generic(){is_face_x(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,-1,0,0) && _const_fm.x > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,-1,0,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,-1,0,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,-1,0,0) < 1e30 ? val(phi,-1,0,0) :
   0.;

 val(ia.x,0,0,0) += _const_alpha.x/(_const_fm.x + 0.)*phif*(val(f,0,0,0) - val(f,-1,0,0))/Delta;
      }}}}end_is_face_x()
#line 89
is_face_y(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,-1,0) && _const_fm.y > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,-1,0) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,-1,0))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,-1,0) < 1e30 ? val(phi,0,-1,0) :
   0.;

 val(ia.y,0,0,0) += _const_alpha.y/(_const_fm.y + 0.)*phif*(val(f,0,0,0) - val(f,0,-1,0))/Delta;
      }}}}end_is_face_y()
#line 89
is_face_z(){
    {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
      if (val(f,0,0,0) != val(f,0,0,-1) && _const_fm.z > 0.) {
#line 101 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
 scalar phi = _attribute[f.i].phi;
 double phif =
   (val(phi,0,0,0) < 1e30 && val(phi,0,0,-1) < 1e30) ?
   (val(phi,0,0,0) + val(phi,0,0,-1))/2. :
   val(phi,0,0,0) < 1e30 ? val(phi,0,0,0) :
   val(phi,0,0,-1) < 1e30 ? val(phi,0,0,-1) :
   0.;

 val(ia.z,0,0,0) += _const_alpha.z/(_const_fm.z + 0.)*phif*(val(f,0,0,0) - val(f,0,0,-1))/Delta;
      }}}}end_is_face_z()}end_foreach_face_generic();}END_FOREACH }
#line 127 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
  {scalar*_i=(scalar*)( list);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){ {
    scalar phi = _attribute[f.i].phi;
    delete (((scalar[]){phi,{-1}}));
    _attribute[f.i].phi.i = 0;
  }}}
  pfree (list,__func__,__FILE__,__LINE__);
}{end_tracing("acceleration_0","/home/esteban/ProgramFile/basilisk/src/iforce.h",133);return 0;}end_tracing("acceleration_0","/home/esteban/ProgramFile/basilisk/src/iforce.h",133);}
#line 16 "/home/esteban/ProgramFile/basilisk/src/tension.h"
#line 1 "curvature.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
#line 68 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
#line 1 "heights.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/heights.h"
#line 29 "/home/esteban/ProgramFile/basilisk/src/heights.h"
static inline double height (double H) {
  return H > 20./2. ? H - 20. : H < -20./2. ? H + 20. : H;
}

static inline int orientation (double H) {
  return fabs(H) > 20./2.;
}
#line 49 "/home/esteban/ProgramFile/basilisk/src/heights.h"
static void half_column (Point point, scalar c, vector h, vector cs, int j)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;






  const int complete = -1;

   {







    double S = val(c,0,0,0), H = S, ci;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.x,0,0,0) == 300.)
 state.s = complete, state.h = 1e30;




      else {
 int s = (val(h.x,0,0,0) + 20./2.)/100.;
 state.h = val(h.x,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/esteban/ProgramFile/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,i*j,0,0) : val(cs.x,(i - 2)*j,0,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      else if (S == ci && trunc(H) != H)
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.x,0,0,0) = 300.;
      else if (S == complete)
 val(h.x,0,0,0) = H;
      else





 val(h.x,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.x,0,0,0) = 1e30;
      else
 val(h.x,0,0,0) = (state.h > 1e10 ? 1e30 : state.h);
    }
  } 
#line 59
{







    double S = val(c,0,0,0), H = S, ci;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.y,0,0,0) == 300.)
 state.s = complete, state.h = 1e30;




      else {
 int s = (val(h.y,0,0,0) + 20./2.)/100.;
 state.h = val(h.y,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/esteban/ProgramFile/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,0,i*j,0) : val(cs.y,0,(i - 2)*j,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      else if (S == ci && trunc(H) != H)
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.y,0,0,0) = 300.;
      else if (S == complete)
 val(h.y,0,0,0) = H;
      else





 val(h.y,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.y,0,0,0) = 1e30;
      else
 val(h.y,0,0,0) = (state.h > 1e10 ? 1e30 : state.h);
    }
  } 
#line 59
{







    double S = val(c,0,0,0), H = S, ci;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.z,0,0,0) == 300.)
 state.s = complete, state.h = 1e30;




      else {
 int s = (val(h.z,0,0,0) + 20./2.)/100.;
 state.h = val(h.z,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/esteban/ProgramFile/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,0,0,i*j) : val(cs.z,0,0,(i - 2)*j);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      else if (S == ci && trunc(H) != H)
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.z,0,0,0) = 300.;
      else if (S == complete)
 val(h.z,0,0,0) = H;
      else





 val(h.z,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.z,0,0,0) = 1e30;
      else
 val(h.z,0,0,0) = (state.h > 1e10 ? 1e30 : state.h);
    }
  }
}
#line 49 "/home/esteban/ProgramFile/basilisk/src/heights.h"
static void _stencil_half_column (Point point, scalar c, vector h, vector cs, int j)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;    






  

   {      







     _stencil_val_o(c,0,0,0);            







    
    
    if (j == 1) {




_stencil_val_o(h.x,0,0,0);{ 
      




{     
 _stencil_val_o(h.x,0,0,0); 
_stencil_val_o(h.x,0,0,0);
    
     
      
#line 92
}}   




         




      





      
      
    
#line 100
}
#line 109 "/home/esteban/ProgramFile/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) { 
_stencil_val_o(c,i*j,0,0); _stencil_val_o(cs.x,(i - 2)*j,0,0); 
           
  
 
       
         
  
 
        







     
   
  
   
        
        
           
        




             
#line 138 "/home/esteban/ProgramFile/basilisk/src/heights.h"
              
              
#line 156 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      
         
    }





    if (j == -1) {







_stencil_val_o(c,0,0,0); _stencil_val_o(c,0,0,0);{
 
{_stencil_val_a(h.x,0,0,0);  }
{
 {_stencil_val_a(h.x,0,0,0);  }





 
{_stencil_val_a(h.x,0,0,0);        }}    
      
#line 183
}







                
              
      
    
#line 184
}
    else {
#line 203
{
 {_stencil_val_a(h.x,0,0,0);  }
 
{_stencil_val_a(h.x,0,0,0);        }}
      
#line 195 "/home/esteban/ProgramFile/basilisk/src/heights.h"
                
   





         
      
    


}
  } 
#line 59
{      







     _stencil_val_o(c,0,0,0);            







    
    
    if (j == 1) {




_stencil_val_o(h.y,0,0,0);{ 
      




{     
 _stencil_val_o(h.y,0,0,0); 
_stencil_val_o(h.y,0,0,0);
    
     
      
#line 92
}}   




         




      





      
      
    
#line 100
}
#line 109 "/home/esteban/ProgramFile/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) { 
_stencil_val_o(c,0,i*j,0); _stencil_val_o(cs.y,0,(i - 2)*j,0); 
           
  
 
       
         
  
 
        







     
   
  
   
        
        
           
        




             
#line 138 "/home/esteban/ProgramFile/basilisk/src/heights.h"
              
              
#line 156 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      
         
    }





    if (j == -1) {







_stencil_val_o(c,0,0,0); _stencil_val_o(c,0,0,0);{
 
{_stencil_val_a(h.y,0,0,0);  }
{
 {_stencil_val_a(h.y,0,0,0);  }





 
{_stencil_val_a(h.y,0,0,0);        }}    
      
#line 183
}







                
              
      
    
#line 184
}
    else {
#line 203
{
 {_stencil_val_a(h.y,0,0,0);  }
 
{_stencil_val_a(h.y,0,0,0);        }}
      
#line 195 "/home/esteban/ProgramFile/basilisk/src/heights.h"
                
   





         
      
    


}
  } 
#line 59
{      







     _stencil_val_o(c,0,0,0);            







    
    
    if (j == 1) {




_stencil_val_o(h.z,0,0,0);{ 
      




{     
 _stencil_val_o(h.z,0,0,0); 
_stencil_val_o(h.z,0,0,0);
    
     
      
#line 92
}}   




         




      





      
      
    
#line 100
}
#line 109 "/home/esteban/ProgramFile/basilisk/src/heights.h"
    for (int i = 1; i <= 4; i++) { 
_stencil_val_o(c,0,0,i*j); _stencil_val_o(cs.z,0,0,(i - 2)*j); 
           
  
 
       
         
  
 
        







     
   
  
   
        
        
           
        




             
#line 138 "/home/esteban/ProgramFile/basilisk/src/heights.h"
              
              
#line 156 "/home/esteban/ProgramFile/basilisk/src/heights.h"
      
         
    }





    if (j == -1) {







_stencil_val_o(c,0,0,0); _stencil_val_o(c,0,0,0);{
 
{_stencil_val_a(h.z,0,0,0);  }
{
 {_stencil_val_a(h.z,0,0,0);  }





 
{_stencil_val_a(h.z,0,0,0);        }}    
      
#line 183
}







                
              
      
    
#line 184
}
    else {
#line 203
{
 {_stencil_val_a(h.z,0,0,0);  }
 
{_stencil_val_a(h.z,0,0,0);        }}
      
#line 195 "/home/esteban/ProgramFile/basilisk/src/heights.h"
                
   





         
      
    


}
  }
}
#line 222 "/home/esteban/ProgramFile/basilisk/src/heights.h"
static void column_propagation (vector h)
{



  foreach_stencil(1,{0},)

    for (int i = -2; i <= 2; i++)
      {
 {_stencil_val(h.x,i,0,0);
_stencil_val(h.x,i,0,0);_stencil_val(h.x,0,0,0);
   { _stencil_val(h.x,i,0,0);_stencil_val_a(h.x,0,0,0);   }      
       
#line 233
}
 
#line 231
{_stencil_val(h.y,0,i,0);
_stencil_val(h.y,0,i,0);_stencil_val(h.y,0,0,0);
   { _stencil_val(h.y,0,i,0);_stencil_val_a(h.y,0,0,0);   }      
       
#line 233
}
 
#line 231
{_stencil_val(h.z,0,0,i);
_stencil_val(h.z,0,0,i);_stencil_val(h.z,0,0,0);
   { _stencil_val(h.z,0,0,i);_stencil_val_a(h.z,0,0,0);   }      
       
#line 233
}}end_foreach_stencil()



   BEGIN_FOREACH{
#line 227
foreach()

    for (int i = -2; i <= 2; i++)
      {
 if (fabs(height(val(h.x,i,0,0))) <= 3.5 &&
     fabs(height(val(h.x,i,0,0)) + i) < fabs(height(val(h.x,0,0,0))))
   val(h.x,0,0,0) = val(h.x,i,0,0) + i;
 
#line 231
if (fabs(height(val(h.y,0,i,0))) <= 3.5 &&
     fabs(height(val(h.y,0,i,0)) + i) < fabs(height(val(h.y,0,0,0))))
   val(h.y,0,0,0) = val(h.y,0,i,0) + i;
 
#line 231
if (fabs(height(val(h.z,0,0,i))) <= 3.5 &&
     fabs(height(val(h.z,0,0,i)) + i) < fabs(height(val(h.z,0,0,0))))
   val(h.z,0,0,0) = val(h.z,0,0,i) + i;}end_foreach();}END_FOREACH 
}
#line 243 "/home/esteban/ProgramFile/basilisk/src/heights.h"
     
void heights (scalar c, vector h)
{tracing("heights","/home/esteban/ProgramFile/basilisk/src/heights.h",244);







  vector  s=new_vector("s");
  
    for (int i = 0; i < nboundary; i++)
      _attribute[s.x.i].boundary[i] = _attribute[c.i].boundary[i];
    
#line 255
for (int i = 0; i < nboundary; i++)
      _attribute[s.y.i].boundary[i] = _attribute[c.i].boundary[i];
    
#line 255
for (int i = 0; i < nboundary; i++)
      _attribute[s.z.i].boundary[i] = _attribute[c.i].boundary[i];






  for (int j = -1; j <= 1; j += 2) {





    foreach_stencil(1,{0},)
      {
        { _stencil_val(c,2*j,0,0);_stencil_val_a(s.x,0,0,0); }
        
#line 271
{ _stencil_val(c,0,2*j,0);_stencil_val_a(s.y,0,0,0); }
        
#line 271
{ _stencil_val(c,0,0,2*j);_stencil_val_a(s.z,0,0,0); }}end_foreach_stencil()





     BEGIN_FOREACH{
#line 269
foreach()
      {
        val(s.x,0,0,0) = val(c,2*j,0,0);
        
#line 271
val(s.y,0,0,0) = val(c,0,2*j,0);
        
#line 271
val(s.z,0,0,0) = val(c,0,0,2*j);}end_foreach();}END_FOREACH 





    foreach_stencil (1,{0},)
      _stencil_half_column (point, c, h, s, j);end_foreach_stencil()





     BEGIN_FOREACH{
#line 277
foreach ()
      half_column (point, c, h, s, j);end_foreach();}END_FOREACH 
  }




  column_propagation (h);delete((scalar*)((vector[]){s,{{-1},{-1},{-1}}}));
end_tracing("heights","/home/esteban/ProgramFile/basilisk/src/heights.h",285);}
#line 459 "/home/esteban/ProgramFile/basilisk/src/heights.h"

#line 69 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
#line 106 "/home/esteban/ProgramFile/basilisk/src/curvature.h"

static double kappa_z (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  int ori = orientation(val(h.z,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h.z,i,j,0) == 1e30 || orientation(val(h.z,i,j,0)) != ori)
 return 1e30;
  double hx = (val(h.z,1,0,0) - val(h.z,-1,0,0))/2.;
  double hy = (val(h.z,0,1,0) - val(h.z,0,-1,0))/2.;







  double filter = 0.2;
  double hxx = (filter*(val(h.z,1,1,0) + val(h.z,-1,1,0) - 2.*val(h.z,0,1,0)) +
  (val(h.z,1,0,0) + val(h.z,-1,0,0) - 2.*val(h.z,0,0,0)) +
  filter*(val(h.z,1,-1,0) + val(h.z,-1,-1,0) - 2.*val(h.z,0,-1,0)))/
    ((1. + 2.*filter)*Delta);
  double hyy = (filter*(val(h.z,1,1,0) + val(h.z,1,-1,0) - 2.*val(h.z,1,0,0)) +
  (val(h.z,0,1,0) + val(h.z,0,-1,0) - 2.*val(h.z,0,0,0)) +
  filter*(val(h.z,-1,1,0) + val(h.z,-1,-1,0) - 2.*val(h.z,-1,0,0)))/
    ((1. + 2.*filter)*Delta);
  double hxy = (val(h.z,1,1,0) + val(h.z,-1,-1,0) - val(h.z,1,-1,0) - val(h.z,-1,1,0))/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}

#line 107
static double kappa_x (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  int ori = orientation(val(h.x,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h.x,0,i,j) == 1e30 || orientation(val(h.x,0,i,j)) != ori)
 return 1e30;
  double hx = (val(h.x,0,1,0) - val(h.x,0,-1,0))/2.;
  double hy = (val(h.x,0,0,1) - val(h.x,0,0,-1))/2.;







  double filter = 0.2;
  double hxx = (filter*(val(h.x,0,1,1) + val(h.x,0,-1,1) - 2.*val(h.x,0,0,1)) +
  (val(h.x,0,1,0) + val(h.x,0,-1,0) - 2.*val(h.x,0,0,0)) +
  filter*(val(h.x,0,1,-1) + val(h.x,0,-1,-1) - 2.*val(h.x,0,0,-1)))/
    ((1. + 2.*filter)*Delta);
  double hyy = (filter*(val(h.x,0,1,1) + val(h.x,0,1,-1) - 2.*val(h.x,0,1,0)) +
  (val(h.x,0,0,1) + val(h.x,0,0,-1) - 2.*val(h.x,0,0,0)) +
  filter*(val(h.x,0,-1,1) + val(h.x,0,-1,-1) - 2.*val(h.x,0,-1,0)))/
    ((1. + 2.*filter)*Delta);
  double hxy = (val(h.x,0,1,1) + val(h.x,0,-1,-1) - val(h.x,0,1,-1) - val(h.x,0,-1,1))/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}

#line 107
static double kappa_y (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  int ori = orientation(val(h.y,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h.y,j,0,i) == 1e30 || orientation(val(h.y,j,0,i)) != ori)
 return 1e30;
  double hx = (val(h.y,0,0,1) - val(h.y,0,0,-1))/2.;
  double hy = (val(h.y,1,0,0) - val(h.y,-1,0,0))/2.;







  double filter = 0.2;
  double hxx = (filter*(val(h.y,1,0,1) + val(h.y,1,0,-1) - 2.*val(h.y,1,0,0)) +
  (val(h.y,0,0,1) + val(h.y,0,0,-1) - 2.*val(h.y,0,0,0)) +
  filter*(val(h.y,-1,0,1) + val(h.y,-1,0,-1) - 2.*val(h.y,-1,0,0)))/
    ((1. + 2.*filter)*Delta);
  double hyy = (filter*(val(h.y,1,0,1) + val(h.y,-1,0,1) - 2.*val(h.y,0,0,1)) +
  (val(h.y,1,0,0) + val(h.y,-1,0,0) - 2.*val(h.y,0,0,0)) +
  filter*(val(h.y,1,0,-1) + val(h.y,-1,0,-1) - 2.*val(h.y,0,0,-1)))/
    ((1. + 2.*filter)*Delta);
  double hxy = (val(h.y,1,0,1) + val(h.y,-1,0,-1) - val(h.y,-1,0,1) - val(h.y,1,0,-1))/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}


static coord normal2_z (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  scalar hz = h.z;
  if (val(hz,0,0,0) == 1e30)
    return (coord){1e30, 1e30, 1e30};
  int ori = orientation(val(hz,0,0,0));
  double a = ori ? -1. : 1.;
  coord n;
  n.z = a;
   {
    if (allocated(-1,0,0) && val(hz,-1,0,0) != 1e30 && orientation(val(hz,-1,0,0)) == ori) {
      if (allocated(1,0,0) && val(hz,1,0,0) != 1e30 && orientation(val(hz,1,0,0)) == ori)
 n.x = a*(val(hz,-1,0,0) - val(hz,1,0,0))/2.;
      else
 n.x = a*(val(hz,-1,0,0) - val(hz,0,0,0));
    }
    else if (allocated(1,0,0) && val(hz,1,0,0) != 1e30 && orientation(val(hz,1,0,0)) == ori)
      n.x = a*(val(hz,0,0,0) - val(hz,1,0,0));
    else
      n.x = 1e30;
  } 
#line 147
{
    if (allocated(0,-1,0) && val(hz,0,-1,0) != 1e30 && orientation(val(hz,0,-1,0)) == ori) {
      if (allocated(0,1,0) && val(hz,0,1,0) != 1e30 && orientation(val(hz,0,1,0)) == ori)
 n.y = a*(val(hz,0,-1,0) - val(hz,0,1,0))/2.;
      else
 n.y = a*(val(hz,0,-1,0) - val(hz,0,0,0));
    }
    else if (allocated(0,1,0) && val(hz,0,1,0) != 1e30 && orientation(val(hz,0,1,0)) == ori)
      n.y = a*(val(hz,0,0,0) - val(hz,0,1,0));
    else
      n.y = 1e30;
  }
  return n;
}

#line 138
static coord normal2_x (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  scalar hz = h.x;
  if (val(hz,0,0,0) == 1e30)
    return (coord){1e30, 1e30, 1e30};
  int ori = orientation(val(hz,0,0,0));
  double a = ori ? -1. : 1.;
  coord n;
  n.x = a;
   {
    if (allocated(0,-1,0) && val(hz,0,-1,0) != 1e30 && orientation(val(hz,0,-1,0)) == ori) {
      if (allocated(0,1,0) && val(hz,0,1,0) != 1e30 && orientation(val(hz,0,1,0)) == ori)
 n.y = a*(val(hz,0,-1,0) - val(hz,0,1,0))/2.;
      else
 n.y = a*(val(hz,0,-1,0) - val(hz,0,0,0));
    }
    else if (allocated(0,1,0) && val(hz,0,1,0) != 1e30 && orientation(val(hz,0,1,0)) == ori)
      n.y = a*(val(hz,0,0,0) - val(hz,0,1,0));
    else
      n.y = 1e30;
  } 
#line 147
{
    if (allocated(0,0,-1) && val(hz,0,0,-1) != 1e30 && orientation(val(hz,0,0,-1)) == ori) {
      if (allocated(0,0,1) && val(hz,0,0,1) != 1e30 && orientation(val(hz,0,0,1)) == ori)
 n.z = a*(val(hz,0,0,-1) - val(hz,0,0,1))/2.;
      else
 n.z = a*(val(hz,0,0,-1) - val(hz,0,0,0));
    }
    else if (allocated(0,0,1) && val(hz,0,0,1) != 1e30 && orientation(val(hz,0,0,1)) == ori)
      n.z = a*(val(hz,0,0,0) - val(hz,0,0,1));
    else
      n.z = 1e30;
  }
  return n;
}

#line 138
static coord normal2_y (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  scalar hz = h.y;
  if (val(hz,0,0,0) == 1e30)
    return (coord){1e30, 1e30, 1e30};
  int ori = orientation(val(hz,0,0,0));
  double a = ori ? -1. : 1.;
  coord n;
  n.y = a;
   {
    if (allocated(0,0,-1) && val(hz,0,0,-1) != 1e30 && orientation(val(hz,0,0,-1)) == ori) {
      if (allocated(0,0,1) && val(hz,0,0,1) != 1e30 && orientation(val(hz,0,0,1)) == ori)
 n.z = a*(val(hz,0,0,-1) - val(hz,0,0,1))/2.;
      else
 n.z = a*(val(hz,0,0,-1) - val(hz,0,0,0));
    }
    else if (allocated(0,0,1) && val(hz,0,0,1) != 1e30 && orientation(val(hz,0,0,1)) == ori)
      n.z = a*(val(hz,0,0,0) - val(hz,0,0,1));
    else
      n.z = 1e30;
  } 
#line 147
{
    if (allocated(-1,0,0) && val(hz,-1,0,0) != 1e30 && orientation(val(hz,-1,0,0)) == ori) {
      if (allocated(1,0,0) && val(hz,1,0,0) != 1e30 && orientation(val(hz,1,0,0)) == ori)
 n.x = a*(val(hz,-1,0,0) - val(hz,1,0,0))/2.;
      else
 n.x = a*(val(hz,-1,0,0) - val(hz,0,0,0));
    }
    else if (allocated(1,0,0) && val(hz,1,0,0) != 1e30 && orientation(val(hz,1,0,0)) == ori)
      n.x = a*(val(hz,0,0,0) - val(hz,1,0,0));
    else
      n.x = 1e30;
  }
  return n;
}


static coord normal_z (Point point, vector h) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  coord n = normal2_z (point, h);
  double nn = fabs(n.x) + fabs(n.y) + fabs(n.z);
  if (nn < 1e30) {
    
      n.x /= nn;
      
#line 168
n.y /= nn;
      
#line 168
n.z /= nn;
    return n;
  }
  return (coord){1e30, 1e30, 1e30};
}

#line 163
static coord normal_x (Point point, vector h) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  coord n = normal2_x (point, h);
  double nn = fabs(n.y) + fabs(n.z) + fabs(n.x);
  if (nn < 1e30) {
    
      n.y /= nn;
      
#line 168
n.z /= nn;
      
#line 168
n.x /= nn;
    return n;
  }
  return (coord){1e30, 1e30, 1e30};
}

#line 163
static coord normal_y (Point point, vector h) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  coord n = normal2_y (point, h);
  double nn = fabs(n.z) + fabs(n.x) + fabs(n.y);
  if (nn < 1e30) {
    
      n.z /= nn;
      
#line 168
n.x /= nn;
      
#line 168
n.y /= nn;
    return n;
  }
  return (coord){1e30, 1e30, 1e30};
}
#line 181 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
static double height_curvature (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;






  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;
  
    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.kappa = kappa_x;
    
#line 195
n.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.kappa = kappa_y;
    
#line 195
n.z.n = val(c,0,0,1) - val(c,0,0,-1), n.z.kappa = kappa_z;
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);




  if (fabs(n.x.n) < fabs(n.y.n))
    do { NormKappa _tmp_ = n.x; n.x = n.y; n.y = _tmp_; } while(false);

  if (fabs(n.x.n) < fabs(n.z.n))
    do { NormKappa _tmp_ = n.x; n.x = n.z; n.z = _tmp_; } while(false);
  if (fabs(n.y.n) < fabs(n.z.n))
    do { NormKappa _tmp_ = n.y; n.y = n.z; n.z = _tmp_; } while(false);





  double kappa = 1e30;
  
    if (kappa == 1e30) {
      kappa = n.x.kappa (point, h);
      if (kappa != 1e30) {
 kappaf = n.x.kappa;
 if (n.x.n < 0.)
   kappa = - kappa;
      }
    }
    
#line 215
if (kappa == 1e30) {
      kappa = n.y.kappa (point, h);
      if (kappa != 1e30) {
 kappaf = n.y.kappa;
 if (n.y.n < 0.)
   kappa = - kappa;
      }
    }
    
#line 215
if (kappa == 1e30) {
      kappa = n.z.kappa (point, h);
      if (kappa != 1e30) {
 kappaf = n.z.kappa;
 if (n.z.n < 0.)
   kappa = - kappa;
      }
    }

  if (kappa != 1e30) {




    if (fabs(kappa) > 1./Delta)
      kappa = sign(kappa)/Delta;
#line 249 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
  }

  return kappa;
}
#line 181 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
static void _stencil_height_curvature (Point point, scalar c,_stencil_undefined * h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;        
      
     
      






  
  
  
    { _stencil_val(c,1,0,0); _stencil_val(c,-1,0,0);     }
    
#line 195
{ _stencil_val(c,0,1,0); _stencil_val(c,0,-1,0);     }
    
#line 195
{ _stencil_val(c,0,0,1); _stencil_val(c,0,0,-1);     }                                                
  
    
    
      




     

     
     





     
     
           
      
          
  
 
      
       
     

  
        




       
#line 249 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
   

  return ;
}






coord height_normal (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;






  typedef struct {
    double n;
    coord (* normal) (Point, vector);
  } NormNormal;
  struct { NormNormal x, y, z; } n;
  
    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.normal = normal_x;
    
#line 273
n.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.normal = normal_y;
    
#line 273
n.z.n = val(c,0,0,1) - val(c,0,0,-1), n.z.normal = normal_z;




  if (fabs(n.x.n) < fabs(n.y.n))
    do { NormNormal _tmp_ = n.x; n.x = n.y; n.y = _tmp_; } while(false);

  if (fabs(n.x.n) < fabs(n.z.n))
    do { NormNormal _tmp_ = n.x; n.x = n.z; n.z = _tmp_; } while(false);
  if (fabs(n.y.n) < fabs(n.z.n))
    do { NormNormal _tmp_ = n.y; n.y = n.z; n.z = _tmp_; } while(false);





  coord normal = {1e30, 1e30, 1e30};
  
    if (normal.x == 1e30)
      normal = n.x.normal (point, h);
    
#line 292
if (normal.y == 1e30)
      normal = n.y.normal (point, h);
    
#line 292
if (normal.z == 1e30)
      normal = n.z.normal (point, h);

  return normal;
}







coord height_normal_z (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  coord nx = normal2_x (point, h);
  coord ny = normal2_y (point, h);
  if (fabs(nx.y) < fabs(ny.x)) {
    normalize (&nx);
    return nx;
  }
  else if (ny.x != 1e30) {
    normalize (&ny);
    return ny;
  }
  return (coord){1e30, 1e30, 1e30};
}

#line 304
coord height_normal_x (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  coord nx = normal2_y (point, h);
  coord ny = normal2_z (point, h);
  if (fabs(nx.z) < fabs(ny.y)) {
    normalize (&nx);
    return nx;
  }
  else if (ny.y != 1e30) {
    normalize (&ny);
    return ny;
  }
  return (coord){1e30, 1e30, 1e30};
}

#line 304
coord height_normal_y (Point point, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  coord nx = normal2_z (point, h);
  coord ny = normal2_x (point, h);
  if (fabs(nx.x) < fabs(ny.z)) {
    normalize (&nx);
    return nx;
  }
  else if (ny.z != 1e30) {
    normalize (&ny);
    return ny;
  }
  return (coord){1e30, 1e30, 1e30};
}
#line 332 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
#line 1 "parabola.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/parabola.h"
#line 1 "utils.h"
#line 2 "/home/esteban/ProgramFile/basilisk/src/parabola.h"






typedef struct {
  coord o;




  double t[3][3];



  double ** M, rhs[6], a[6];


} ParabolaFit;

static void parabola_fit_init (ParabolaFit * p, coord o, coord m)
{
  
    p->o.x = o.x;
    
#line 26
p->o.y = o.y;
    
#line 26
p->o.z = o.z;






  double max;
  coord nx = {0., 0., 0.}, ny, nz;
  int d = 0;

  
    nz.x = m.x;
    
#line 38
nz.y = m.y;
    
#line 38
nz.z = m.z;
  normalize (&nz);
  max = sq(nz.x);

  if (sq(nz.y) > max) { max = sq(nz.y); d = 1; }
  if (sq(nz.z) > max) d = 2;
  switch (d) {
  case 0: nx.x = - nz.z/nz.x; nx.z = 1.0; break;
  case 1: nx.y = - nz.z/nz.y; nx.z = 1.0; break;
  case 2: nx.z = - nz.x/nz.z; nx.x = 1.0; break;
  }
  normalize (&nx);


  
    ny.x = nz.y*nx.z - nz.z*nx.y;
    
#line 53
ny.y = nz.z*nx.x - nz.x*nx.z;
    
#line 53
ny.z = nz.x*nx.y - nz.y*nx.x;


  p->t[0][0] = nx.x; p->t[0][1] = nx.y; p->t[0][2] = nx.z;
  p->t[1][0] = ny.x; p->t[1][1] = ny.y; p->t[1][2] = ny.z;
  p->t[2][0] = nz.x; p->t[2][1] = nz.y; p->t[2][2] = nz.z;



  int n = 6;


  p->M = (double **) matrix_new (n, n, sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      p->M[i][j] = 0.;
    p->rhs[i] = 0.;
  }
}

static void parabola_fit_add (ParabolaFit * p, coord m, double w)
{
#line 85 "/home/esteban/ProgramFile/basilisk/src/parabola.h"
  double x1 = m.x - p->o.x, y1 = m.y - p->o.y, z1 = m.z - p->o.z;
  double x = p->t[0][0]*x1 + p->t[0][1]*y1 + p->t[0][2]*z1;
  double y = p->t[1][0]*x1 + p->t[1][1]*y1 + p->t[1][2]*z1;
  double z = p->t[2][0]*x1 + p->t[2][1]*y1 + p->t[2][2]*z1;
#line 98 "/home/esteban/ProgramFile/basilisk/src/parabola.h"
  double x2 = x*x, x3 = x2*x, x4 = x3*x;
  double y2 = y*y, y3 = y2*y, y4 = y3*y;
  p->M[0][0] += w*x4; p->M[1][1] += w*y4; p->M[2][2] += w*x2*y2;
  p->M[3][3] += w*x2; p->M[4][4] += w*y2; p->M[5][5] += w;
  p->M[0][2] += w*x3*y; p->M[0][3] += w*x3; p->M[0][4] += w*x2*y;
  p->M[1][2] += w*x*y3; p->M[1][3] += w*x*y2; p->M[1][4] += w*y3;
  p->M[2][5] += w*x*y;
  p->M[3][5] += w*x;
  p->M[4][5] += w*y;
  p->rhs[0] += w*x2*z; p->rhs[1] += w*y2*z; p->rhs[2] += w*x*y*z;
  p->rhs[3] += w*x*z; p->rhs[4] += w*y*z; p->rhs[5] += w*z;


}

static double parabola_fit_solve (ParabolaFit * p)
{
#line 139 "/home/esteban/ProgramFile/basilisk/src/parabola.h"
  p->M[0][1] = p->M[2][2]; p->M[0][5] = p->M[3][3];
  p->M[1][5] = p->M[4][4];
  p->M[2][3] = p->M[0][4]; p->M[2][4] = p->M[1][3];
  p->M[3][4] = p->M[2][5];
  for (int i = 1; i < 6; i++)
    for (int j = 0; j < i; j++)
      p->M[i][j] = p->M[j][i];
  double pivmin = matrix_inverse (p->M, 6, 1e-10);
  if (pivmin)
    for (int i = 0; i < 6; i++) {
      p->a[i] = 0.;
      for (int j = 0; j < 6; j++)
 p->a[i] += p->M[i][j]*p->rhs[j];
    }
  else
    for (int i = 0; i < 6; i++)
      p->a[i] = 0.;


  matrix_free (p->M);
  return pivmin;
}

static double parabola_fit_curvature (ParabolaFit * p,
          double kappamax, double * kmax)
{
  double kappa;
#line 176 "/home/esteban/ProgramFile/basilisk/src/parabola.h"
  double hxx = 2.*p->a[0], hyy = 2.*p->a[1], hxy = p->a[2];
  double hx = p->a[3], hy = p->a[4];

  double dnm = 1. + sq(hx) + sq(hy);
  kappa = - (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)
    /sqrt (dnm*dnm*dnm);
  if (kmax) {
    double kg = (hxx*hyy - hxy*hxy)/(dnm*dnm);
    double a = kappa*kappa/4. - kg;
    *kmax = fabs (kappa/2.);
    if (a >= 0.)
      *kmax += sqrt (a);
  }

  if (fabs (kappa) > kappamax) {
    if (kmax)
      *kmax = kappamax;
    return kappa > 0. ? kappamax : - kappamax;
  }
  return kappa;
}
#line 333 "/home/esteban/ProgramFile/basilisk/src/curvature.h"






static int independents (coord * p, int n)
{
  if (n < 2)
    return n;
  int ni = 1;
  for (int j = 1; j < n; j++) {
    bool depends = false;
    for (int i = 0; i < j && !depends; i++) {
      double d2 = 0.;
      
 d2 += sq(p[i].x - p[j].x);
 
#line 349
d2 += sq(p[i].y - p[j].y);
 
#line 349
d2 += sq(p[i].z - p[j].z);
      depends = (d2 < sq(0.5));
    }
    ni += !depends;
  }
  return ni;
}






static double height_curvature_fit (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;





  coord ip[3 == 2 ? 6 : 27];
  int n = 0;




   {





    int n1 = 0, n2 = 0;






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.z,i,j,0) != 1e30) {
   if (orientation(val(h.z,i,j,0))) n1++; else n2++;
 }

    int ori = (n1 > n2);
#line 406 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.z,i,j,0) != 1e30 && orientation(val(h.z,i,j,0)) == ori)
   ip[n].x = i, ip[n].y = j, ip[n++].z = height(val(h.z,i,j,0));

  } 
#line 375
{





    int n1 = 0, n2 = 0;






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.x,0,i,j) != 1e30) {
   if (orientation(val(h.x,0,i,j))) n1++; else n2++;
 }

    int ori = (n1 > n2);
#line 406 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.x,0,i,j) != 1e30 && orientation(val(h.x,0,i,j)) == ori)
   ip[n].y = i, ip[n].z = j, ip[n++].x = height(val(h.x,0,i,j));

  } 
#line 375
{





    int n1 = 0, n2 = 0;






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.y,j,0,i) != 1e30) {
   if (orientation(val(h.y,j,0,i))) n1++; else n2++;
 }

    int ori = (n1 > n2);
#line 406 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.y,j,0,i) != 1e30 && orientation(val(h.y,j,0,i)) == ori)
   ip[n].z = i, ip[n].x = j, ip[n++].y = height(val(h.y,j,0,i));

  }





  if (independents (ip, n) < (3 == 2 ? 3 : 9))
    return 1e30;





  coord m = mycs (point, c), fc;
  double alpha = plane_alpha (val(c,0,0,0), m);
  double area = plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);




  parabola_fit_add (&fit, fc, area*100.);






  for (int i = 0; i < n; i++)
    parabola_fit_add (&fit, ip[i], 1.);
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}







#line 362
static void _stencil_height_curvature_fit (Point point, scalar c, vector h)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;          





  
  




   {      





    






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 {_stencil_val(h.z,i,j,0); {
_stencil_val(h.z,i,j,0);  
     
 
#line 392
}   }     

    
#line 406 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 {_stencil_val(h.z,i,j,0);_stencil_val(h.z,i,j,0);
   {_stencil_val(h.z,i,j,0);        }       }

  } 
#line 375
{      





    






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 {_stencil_val(h.x,0,i,j); {
_stencil_val(h.x,0,i,j);  
     
 
#line 392
}   }     

    
#line 406 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 {_stencil_val(h.x,0,i,j);_stencil_val(h.x,0,i,j);
   {_stencil_val(h.x,0,i,j);        }       }

  } 
#line 375
{      





    






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 {_stencil_val(h.y,j,0,i); {
_stencil_val(h.y,j,0,i);  
     
 
#line 392
}   }     

    
#line 406 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 {_stencil_val(h.y,j,0,i);_stencil_val(h.y,j,0,i);
   {_stencil_val(h.y,j,0,i);        }       }

  }    
    





             





   _stencil_mycs (point, c);     
  _stencil_val(c,0,0,0);             
  
                 
  




  






     
    
  
  



  return ;
}






static double centroids_curvature_fit (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;





  coord m = mycs (point, c), fc;
  double alpha = plane_alpha (val(c,0,0,0), m);
  plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);





  coord r = {x,y,z};
  {foreach_neighbor(1)
    if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
      coord m = mycs (point, c), fc;
      double alpha = plane_alpha (val(c,0,0,0), m);
      double area = plane_area_center (m, alpha, &fc);
      coord rn = {x,y,z};
      
 fc.x += (rn.x - r.x)/Delta;
 
#line 480
fc.y += (rn.y - r.y)/Delta;
 
#line 480
fc.z += (rn.z - r.z)/Delta;
      parabola_fit_add (&fit, fc, area);
    }end_foreach_neighbor()}
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}







#line 455
static void _stencil_centroids_curvature_fit (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;   





   _stencil_mycs (point, c);     
  _stencil_val(c,0,0,0);       
  
     
  





  
  {foreach_neighbor(1)
    {_stencil_val(c,0,0,0); _stencil_val(c,0,0,0); {   
       _stencil_mycs (point, c);     
      _stencil_val(c,0,0,0);         
      
      
       
       
      
    }      }end_foreach_neighbor()}       
  
  



  return ;
}
#line 504 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
static inline bool interfacial (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if (val(c,0,0,0) >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      {
 if (val(c,i,0,0) <= 0.)
   return true;
 
#line 509
if (val(c,0,i,0) <= 0.)
   return true;
 
#line 509
if (val(c,0,0,i) <= 0.)
   return true;}
  }
  else if (val(c,0,0,0) <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      {
 if (val(c,i,0,0) >= 1.)
   return true;
 
#line 515
if (val(c,0,i,0) >= 1.)
   return true;
 
#line 515
if (val(c,0,0,i) >= 1.)
   return true;}
  }
  else
    return true;
  return false;
}
#line 504 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
static void _stencil_interfacial (Point point, scalar c)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
_stencil_val(c,0,0,0);{ {
    for (int i = -1; i <= 1; i += 2)
      {
 {_stencil_val(c,i,0,0); 
      }
 
#line 509
{_stencil_val(c,0,i,0); 
      }
 
#line 509
{_stencil_val(c,0,0,i); 
      }}
  } 
{_stencil_val(c,0,0,0);{ {
    for (int i = -1; i <= 1; i += 2)
      {
 {_stencil_val(c,i,0,0); 
      }
 
#line 515
{_stencil_val(c,0,i,0); 
      }
 
#line 515
{_stencil_val(c,0,0,i); 
      }}
  } 
    
}   
  
#line 519
}}
     
  
  
#line 520
return ;
}
#line 533 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
typedef struct {
  int h;
  int f;
  int a;
  int c;
} cstats;

     
cstats curvature (scalar c, scalar kappa,
    double sigma, bool add)
{tracing("curvature","/home/esteban/ProgramFile/basilisk/src/curvature.h",541);
  int sh = 0, f = 0, sa = 0, sc = 0;
#line 557 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
  vector ch = _attribute[c.i].height,   h=(ch).x.i>0?(ch):new_vector("h");
  if (!ch.x.i)
    heights (c, h);





  scalar  k=new_scalar("k");
  scalar_clone (k, kappa);

  foreach_stencil(1,{0},) {




_stencil_interfacial (point, c);{
      {_stencil_val_a(k,0,0,0);  } 





{ _stencil_height_curvature (point, c, NULL);_stencil_val_a(k,0,0,0);{
       
{ _stencil_height_curvature_fit (point, c, h);_stencil_val_a(k,0,0,0);
          }}    
    
#line 583
}}




     





    
  
#line 584
}end_foreach_stencil()

  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(+:f)reduction(+:sh)){
#line 568
foreach() {




    if (!interfacial (point, c))
      val(k,0,0,0) = 1e30;





    else if ((val(k,0,0,0) = height_curvature (point, c, h)) != 1e30)
      sh++;
    else if ((val(k,0,0,0) = height_curvature_fit (point, c, h)) != 1e30)
      f++;
  }end_foreach();mpi_all_reduce_array(&f,int,MPI_SUM,1);mpi_all_reduce_array(&sh,int,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 

  
#line 586
foreach_stencil (1,{0},) { 





    
_stencil_val(k,0,0,0);{
      { _stencil_val(k,0,0,0); } 
{_stencil_interfacial (point, c);{ {      





      
      {foreach_neighbor(1)
 {_stencil_val(k,0,0,0);
   { _stencil_val(k,0,0,0);  }   }end_foreach_neighbor()}




 


{ _stencil_centroids_curvature_fit (point, c);  }
         
    
      
    
#line 613
}
        
} 
    
#line 615
}}




{
      {_stencil_val_a(kappa,0,0,0);  } 
if (add)
      {_stencil_val_r(kappa,0,0,0);  }
    else
      {_stencil_val_a(kappa,0,0,0);  }}
       
    




       
    
  
#line 626
}end_foreach_stencil()

  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(+:sc)reduction(+:sa)){
#line 586
foreach () {





    double kf;
    if (val(k,0,0,0) < 1e30)
      kf = val(k,0,0,0);
    else if (interfacial (point, c)) {





      double sk = 0., a = 0.;
      {foreach_neighbor(1)
 if (val(k,0,0,0) < 1e30)
   sk += val(k,0,0,0), a++;end_foreach_neighbor()}
      if (a > 0.)
 kf = sk/a, sa++;
      else




 kf = centroids_curvature_fit (point, c), sc++;
    }
    else
      kf = 1e30;




    if (kf == 1e30)
      val(kappa,0,0,0) = 1e30;
    else if (add)
      val(kappa,0,0,0) += sigma*kf;
    else
      val(kappa,0,0,0) = sigma*kf;
  }end_foreach();mpi_all_reduce_array(&sc,int,MPI_SUM,1);mpi_all_reduce_array(&sa,int,MPI_SUM,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 643 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
  { cstats _ret= (cstats){sh, f, sa, sc};delete((scalar*)((scalar[]){k,{-1}}));if((ch).x.i<=0)delete((scalar*)((vector[]){h,{{-1},{-1},{-1}}}));{end_tracing("curvature","/home/esteban/ProgramFile/basilisk/src/curvature.h",643);return _ret;}}delete((scalar*)((scalar[]){k,{-1}}));
end_tracing("curvature","/home/esteban/ProgramFile/basilisk/src/curvature.h",644);}
#line 665 "/home/esteban/ProgramFile/basilisk/src/curvature.h"

static double pos_x (Point point, vector h, coord * G, coord * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if (fabs(height(val(h.x,0,0,0))) > 1.)
    return 1e30;
  coord o = {x, y, z};
  o.x += height(val(h.x,0,0,0))*Delta;
  double pos = 0.;
  
    pos += (o.x - Z->x)*G->x;
    
#line 674
pos += (o.y - Z->y)*G->y;
    
#line 674
pos += (o.z - Z->z)*G->z;
  return pos;
}

#line 666
static double pos_y (Point point, vector h, coord * G, coord * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if (fabs(height(val(h.y,0,0,0))) > 1.)
    return 1e30;
  coord o = {x, y, z};
  o.y += height(val(h.y,0,0,0))*Delta;
  double pos = 0.;
  
    pos += (o.y - Z->y)*G->y;
    
#line 674
pos += (o.z - Z->z)*G->z;
    
#line 674
pos += (o.x - Z->x)*G->x;
  return pos;
}

#line 666
static double pos_z (Point point, vector h, coord * G, coord * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if (fabs(height(val(h.z,0,0,0))) > 1.)
    return 1e30;
  coord o = {x, y, z};
  o.z += height(val(h.z,0,0,0))*Delta;
  double pos = 0.;
  
    pos += (o.z - Z->z)*G->z;
    
#line 674
pos += (o.x - Z->x)*G->x;
    
#line 674
pos += (o.y - Z->y)*G->y;
  return pos;
}







static double height_position (Point point, scalar f, vector h,
          coord * G, coord * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;






  typedef struct {
    double n;
    double (* pos) (Point, vector, coord *, coord *);
  } NormPos;
  struct { NormPos x, y, z; } n;
  
    n.x.n = val(f,1,0,0) - val(f,-1,0,0), n.x.pos = pos_x;
    
#line 699
n.y.n = val(f,0,1,0) - val(f,0,-1,0), n.y.pos = pos_y;
    
#line 699
n.z.n = val(f,0,0,1) - val(f,0,0,-1), n.z.pos = pos_z;




  if (fabs(n.x.n) < fabs(n.y.n))
    do { NormPos _tmp_ = n.x; n.x = n.y; n.y = _tmp_; } while(false);

  if (fabs(n.x.n) < fabs(n.z.n))
    do { NormPos _tmp_ = n.x; n.x = n.z; n.z = _tmp_; } while(false);
  if (fabs(n.y.n) < fabs(n.z.n))
    do { NormPos _tmp_ = n.y; n.y = n.z; n.z = _tmp_; } while(false);





  double pos = 1e30;
  
    if (pos == 1e30)
      pos = n.x.pos (point, h, G, Z);
    
#line 718
if (pos == 1e30)
      pos = n.y.pos (point, h, G, Z);
    
#line 718
if (pos == 1e30)
      pos = n.z.pos (point, h, G, Z);

  return pos;
}








#line 684
static void _stencil_height_position (Point point, scalar f,_stencil_undefined * h,
_stencil_undefined * G,_stencil_undefined * Z)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;        
          
     
      






  
  
  
    { _stencil_val(f,1,0,0); _stencil_val(f,-1,0,0);     }
    
#line 699
{ _stencil_val(f,0,1,0); _stencil_val(f,0,-1,0);     }
    
#line 699
{ _stencil_val(f,0,0,1); _stencil_val(f,0,0,-1);     }                                          
    
    
    




     

     
     





  
     
          
      

  return ;
}
#line 735 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
void position (scalar f, scalar pos,
        coord G, coord Z, bool add)
{
#line 749 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
  vector fh = _attribute[f.i].height,   h=(fh).x.i>0?(fh):new_vector("h");
  if (!fh.x.i)
    heights (f, h);
  foreach_stencil(1,{0},) {
_stencil_interfacial (point, f);{ {  
       _stencil_height_position (point, f, NULL,NULL ,NULL ); 
{      





  _stencil_mycs (point, f);     
 _stencil_val(f,0,0,0);    
 
  
 
         
      }
         
      
#line 768
if (add)
 {_stencil_val_r(pos,0,0,0);  }
      else
 {_stencil_val_a(pos,0,0,0);  }
    }
      
{_stencil_val_a(pos,0,0,0);  }}
     
    
  
#line 775
}end_foreach_stencil()
   BEGIN_FOREACH{
#line 752
foreach() {
    if (interfacial (point, f)) {
      double hp = height_position (point, f, h, &G, &Z);
      if (hp == 1e30) {





 coord n = mycs (point, f), o = {x,y,z}, c;
 double alpha = plane_alpha (val(f,0,0,0), n);
 plane_area_center (n, alpha, &c);
 hp = 0.;
 
   hp += (o.x + Delta*c.x - Z.x)*G.x;
   
#line 766
hp += (o.y + Delta*c.y - Z.y)*G.y;
   
#line 766
hp += (o.z + Delta*c.z - Z.z)*G.z;
      }
      if (add)
 val(pos,0,0,0) += hp;
      else
 val(pos,0,0,0) = hp;
    }
    else
      val(pos,0,0,0) = 1e30;
  }end_foreach();}END_FOREACH if((fh).x.i<=0)delete((scalar*)((vector[]){h,{{-1},{-1},{-1}}}));
#line 790 "/home/esteban/ProgramFile/basilisk/src/curvature.h"
}
#line 36 "/home/esteban/ProgramFile/basilisk/src/tension.h"
static int stability_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 17 "/home/esteban/ProgramFile/basilisk/src/tension.h"





#line 36 "/home/esteban/ProgramFile/basilisk/src/tension.h"
      static int stability_1(const int i,const double t,Event *_ev){tracing("stability_1","/home/esteban/ProgramFile/basilisk/src/tension.h",36);
{





  double amin = 1e30, amax = -1e30, dmin = 1e30;
  if(!is_constant(fm.x) && !is_constant(alpha.x)){
  
#line 44
foreach_face_stencil (1,{0},){_stencil_is_face_x(){
    {_stencil_val(fm.x,0,0,0); {
_stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); { _stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); }
_stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); { _stencil_val(alpha.x,0,0,0);_stencil_val(fm.x,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_x()
#line 44
_stencil_is_face_y(){
    {_stencil_val(fm.y,0,0,0); {
_stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); { _stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); }
_stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); { _stencil_val(alpha.y,0,0,0);_stencil_val(fm.y,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_y()
#line 44
_stencil_is_face_z(){
    {_stencil_val(fm.z,0,0,0); {
_stencil_val(alpha.z,0,0,0);_stencil_val(fm.z,0,0,0); { _stencil_val(alpha.z,0,0,0);_stencil_val(fm.z,0,0,0); }
_stencil_val(alpha.z,0,0,0);_stencil_val(fm.z,0,0,0); { _stencil_val(alpha.z,0,0,0);_stencil_val(fm.z,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_z()}end_foreach_face_stencil()
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (val(fm.x,0,0,0) > 0.) {
      if (val(alpha.x,0,0,0)/val(fm.x,0,0,0) > amax) amax = val(alpha.x,0,0,0)/val(fm.x,0,0,0);
      if (val(alpha.x,0,0,0)/val(fm.x,0,0,0) < amin) amin = val(alpha.x,0,0,0)/val(fm.x,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (val(fm.y,0,0,0) > 0.) {
      if (val(alpha.y,0,0,0)/val(fm.y,0,0,0) > amax) amax = val(alpha.y,0,0,0)/val(fm.y,0,0,0);
      if (val(alpha.y,0,0,0)/val(fm.y,0,0,0) < amin) amin = val(alpha.y,0,0,0)/val(fm.y,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()
#line 44
is_face_z(){
    if (val(fm.z,0,0,0) > 0.) {
      if (val(alpha.z,0,0,0)/val(fm.z,0,0,0) > amax) amax = val(alpha.z,0,0,0)/val(fm.z,0,0,0);
      if (val(alpha.z,0,0,0)/val(fm.z,0,0,0) < amin) amin = val(alpha.z,0,0,0)/val(fm.z,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_z()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 49
}else if(is_constant(fm.x) && !is_constant(alpha.x)){_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);
  
#line 44
foreach_face_stencil (1,{0},){_stencil_is_face_x(){
    {; {
_stencil_val(alpha.x,0,0,0);; { _stencil_val(alpha.x,0,0,0);; }
_stencil_val(alpha.x,0,0,0);; { _stencil_val(alpha.x,0,0,0);; }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_x()
#line 44
_stencil_is_face_y(){
    {; {
_stencil_val(alpha.y,0,0,0);; { _stencil_val(alpha.y,0,0,0);; }
_stencil_val(alpha.y,0,0,0);; { _stencil_val(alpha.y,0,0,0);; }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_y()
#line 44
_stencil_is_face_z(){
    {; {
_stencil_val(alpha.z,0,0,0);; { _stencil_val(alpha.z,0,0,0);; }
_stencil_val(alpha.z,0,0,0);; { _stencil_val(alpha.z,0,0,0);; }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_z()}end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (_const_fm.x > 0.) {
      if (val(alpha.x,0,0,0)/_const_fm.x > amax) amax = val(alpha.x,0,0,0)/_const_fm.x;
      if (val(alpha.x,0,0,0)/_const_fm.x < amin) amin = val(alpha.x,0,0,0)/_const_fm.x;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (_const_fm.y > 0.) {
      if (val(alpha.y,0,0,0)/_const_fm.y > amax) amax = val(alpha.y,0,0,0)/_const_fm.y;
      if (val(alpha.y,0,0,0)/_const_fm.y < amin) amin = val(alpha.y,0,0,0)/_const_fm.y;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()
#line 44
is_face_z(){
    if (_const_fm.z > 0.) {
      if (val(alpha.z,0,0,0)/_const_fm.z > amax) amax = val(alpha.z,0,0,0)/_const_fm.z;
      if (val(alpha.z,0,0,0)/_const_fm.z < amin) amin = val(alpha.z,0,0,0)/_const_fm.z;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_z()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 49
}else if(!is_constant(fm.x) && is_constant(alpha.x)){_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 44
foreach_face_stencil (1,{0},){_stencil_is_face_x(){
    {_stencil_val(fm.x,0,0,0); {
;_stencil_val(fm.x,0,0,0); {;_stencil_val(fm.x,0,0,0); }
;_stencil_val(fm.x,0,0,0); {;_stencil_val(fm.x,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_x()
#line 44
_stencil_is_face_y(){
    {_stencil_val(fm.y,0,0,0); {
;_stencil_val(fm.y,0,0,0); {;_stencil_val(fm.y,0,0,0); }
;_stencil_val(fm.y,0,0,0); {;_stencil_val(fm.y,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_y()
#line 44
_stencil_is_face_z(){
    {_stencil_val(fm.z,0,0,0); {
;_stencil_val(fm.z,0,0,0); {;_stencil_val(fm.z,0,0,0); }
;_stencil_val(fm.z,0,0,0); {;_stencil_val(fm.z,0,0,0); }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_z()}end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (val(fm.x,0,0,0) > 0.) {
      if (_const_alpha.x/val(fm.x,0,0,0) > amax) amax = _const_alpha.x/val(fm.x,0,0,0);
      if (_const_alpha.x/val(fm.x,0,0,0) < amin) amin = _const_alpha.x/val(fm.x,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (val(fm.y,0,0,0) > 0.) {
      if (_const_alpha.y/val(fm.y,0,0,0) > amax) amax = _const_alpha.y/val(fm.y,0,0,0);
      if (_const_alpha.y/val(fm.y,0,0,0) < amin) amin = _const_alpha.y/val(fm.y,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()
#line 44
is_face_z(){
    if (val(fm.z,0,0,0) > 0.) {
      if (_const_alpha.z/val(fm.z,0,0,0) > amax) amax = _const_alpha.z/val(fm.z,0,0,0);
      if (_const_alpha.z/val(fm.z,0,0,0) < amin) amin = _const_alpha.z/val(fm.z,0,0,0);
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_z()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 49
}else {_coord _const_fm={_constant[fm.x.i-_NVARMAX],_constant[fm.y.i-_NVARMAX],_constant[fm.z.i-_NVARMAX]};NOT_UNUSED(_const_fm);_coord _const_alpha={_constant[alpha.x.i-_NVARMAX],_constant[alpha.y.i-_NVARMAX],_constant[alpha.z.i-_NVARMAX]};NOT_UNUSED(_const_alpha);
  
#line 44
foreach_face_stencil (1,{0},){_stencil_is_face_x(){
    {; {
;; {;; }
;; {;; }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_x()
#line 44
_stencil_is_face_y(){
    {; {
;; {;; }
;; {;; }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_y()
#line 44
_stencil_is_face_z(){
    {; {
;; {;; }
;; {;; }   
         
         
         
    
#line 49
}   }}end__stencil_is_face_z()}end_foreach_face_stencil()
  
#undef OMP_PARALLEL
#define OMP_PARALLEL()
BEGIN_FOREACH OMP(omp parallel  reduction(min:dmin) reduction(max:amax)reduction(min:amin)){
#line 44
foreach_face_generic (){is_face_x(){
    if (_const_fm.x > 0.) {
      if (_const_alpha.x/_const_fm.x > amax) amax = _const_alpha.x/_const_fm.x;
      if (_const_alpha.x/_const_fm.x < amin) amin = _const_alpha.x/_const_fm.x;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_x()
#line 44
is_face_y(){
    if (_const_fm.y > 0.) {
      if (_const_alpha.y/_const_fm.y > amax) amax = _const_alpha.y/_const_fm.y;
      if (_const_alpha.y/_const_fm.y < amin) amin = _const_alpha.y/_const_fm.y;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_y()
#line 44
is_face_z(){
    if (_const_fm.z > 0.) {
      if (_const_alpha.z/_const_fm.z > amax) amax = _const_alpha.z/_const_fm.z;
      if (_const_alpha.z/_const_fm.z < amin) amin = _const_alpha.z/_const_fm.z;
      if (Delta < dmin) dmin = Delta;
    }}end_is_face_z()}end_foreach_face_generic();mpi_all_reduce_array(&dmin,double,MPI_MIN,1);mpi_all_reduce_array(&amax,double,MPI_MAX,1);mpi_all_reduce_array(&amin,double,MPI_MIN,1);
#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}END_FOREACH 
#line 49
}
  double rhom = (1./amin + 1./amax)/2.;





  double sigma = 0.;
  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar c=*_i;(&c)->i>=0;c=*++_i){
    sigma += _attribute[c.i].sigma;}}
  if (sigma) {
    double dt = sqrt (rhom*cube(dmin)/(3.14159265358979*sigma));
    if (dt < dtmax)
      dtmax = dt;
  }
}{end_tracing("stability_1","/home/esteban/ProgramFile/basilisk/src/tension.h",64);return 0;}end_tracing("stability_1","/home/esteban/ProgramFile/basilisk/src/tension.h",64);}







static int acceleration_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}








#line 72
      static int acceleration_1(const int i,const double t,Event *_ev){tracing("acceleration_1","/home/esteban/ProgramFile/basilisk/src/tension.h",72);
{




  {scalar*_i=(scalar*)( interfaces);if(_i)for(scalar f=*_i;(&f)->i>=0;f=*++_i){
    if (_attribute[f.i].sigma) {





      scalar phi = _attribute[f.i].phi;
      if (phi.i)
 curvature (f, phi, _attribute[f.i].sigma, true);
      else {
 phi = new_scalar("phi");
 curvature (f, phi, _attribute[f.i].sigma, false);
 _attribute[f.i].phi = phi;
      }
    }}}
}{end_tracing("acceleration_1","/home/esteban/ProgramFile/basilisk/src/tension.h",94);return 0;}end_tracing("acceleration_1","/home/esteban/ProgramFile/basilisk/src/tension.h",94);}
#line 15 "main3D.c"
#line 1 "reduced.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/reduced.h"
#line 19 "/home/esteban/ProgramFile/basilisk/src/reduced.h"
coord G = {0.,0.,0.}, Z = {0.,0.,0.};
#line 36
static int acceleration_2_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}
#line 36 "/home/esteban/ProgramFile/basilisk/src/reduced.h"
      static int acceleration_2(const int i,const double t,Event *_ev){tracing("acceleration_2","/home/esteban/ProgramFile/basilisk/src/reduced.h",36);
{
  scalar phi = _attribute[f.i].phi;
  coord G1;
  
    G1.x = (rho2 - rho1)*G.x;
    
#line 41
G1.y = (rho2 - rho1)*G.y;
    
#line 41
G1.z = (rho2 - rho1)*G.z;

  if (phi.i)
    position (f, phi, G1, Z, true);
  else {
    phi = new_scalar("phi");
    position (f, phi, G1, Z, false);
    _attribute[f.i].phi = phi;
  }
}{end_tracing("acceleration_2","/home/esteban/ProgramFile/basilisk/src/reduced.h",50);return 0;}end_tracing("acceleration_2","/home/esteban/ProgramFile/basilisk/src/reduced.h",50);}
#line 37 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h"
static int defaults_4_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i = 0)!=0;*ip=i;*tp=t;return ret;}
#line 16 "main3D.c"

#line 1 "navier-stokes/conserving.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h"
#line 37 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h"
      static int defaults_4(const int i,const double t,Event *_ev){tracing("defaults_4","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",37);
{
  stokes = true;
#line 66 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h"
}{end_tracing("defaults_4","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",66);return 0;}end_tracing("defaults_4","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",66);}





static int stability_2_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}






#line 72
      static int stability_2(const int i,const double t,Event *_ev){tracing("stability_2","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",72);
{dtmax = timestep (uf, dtmax);
  
#line 73
}{end_tracing("stability_2","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",73);return 0;}end_tracing("stability_2","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",73);}








static double boundary_q1_x (Point neighbor, Point point, scalar q1, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return clamp(val(f,0,0,0),0.,1.)*rho1*val(u.x,0,0,0);
}

#line 82
static double boundary_q1_y (Point neighbor, Point point, scalar q1, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return clamp(val(f,0,0,0),0.,1.)*rho1*val(u.y,0,0,0);
}

#line 82
static double boundary_q1_z (Point neighbor, Point point, scalar q1, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return clamp(val(f,0,0,0),0.,1.)*rho1*val(u.z,0,0,0);
}


static double boundary_q2_x (Point neighbor, Point point, scalar q2, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.x,0,0,0);
}

#line 88
static double boundary_q2_y (Point neighbor, Point point, scalar q2, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.y,0,0,0);
}

#line 88
static double boundary_q2_z (Point neighbor, Point point, scalar q2, bool * data)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return (1. - clamp(val(f,0,0,0),0.,1.))*rho2*val(u.z,0,0,0);
}
#line 115 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h"
static scalar * interfaces1 = NULL;

static int vof_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}


#line 117
      static int vof_1(const int i,const double t,Event *_ev){tracing("vof_1","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",117); {







  vector  q1=new_vector("q1"),  q2=new_vector("q2");
  {scalar*_i=(scalar*)(((vector[]) {q1,q2,{{-1},{-1},{-1}}}));if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
    _attribute[s.i].depends = list_add (_attribute[s.i].depends, f);
    
      _attribute[s.i].v.x.i = -1;
      
#line 129
_attribute[s.i].v.y.i = -1;
      
#line 129
_attribute[s.i].v.z.i = -1;
  }}}
  for (int i = 0; i < nboundary; i++)
    { {
      _attribute[q1.x.i].boundary[i] = boundary_q1_x;
      _attribute[q2.x.i].boundary[i] = boundary_q2_x;
    } 
#line 132
{
      _attribute[q1.y.i].boundary[i] = boundary_q1_y;
      _attribute[q2.y.i].boundary[i] = boundary_q2_y;
    } 
#line 132
{
      _attribute[q1.z.i].boundary[i] = boundary_q1_z;
      _attribute[q2.z.i].boundary[i] = boundary_q2_z;
    }}
#line 147 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h"
  foreach_stencil(1,{0},)
    { {   
      _stencil_val(f,0,0,0);
_stencil_val(u.x,0,0,0);
      
#line 150
_stencil_val_a(q1.x,0,0,0);
_stencil_val(u.x,0,0,0);  
      
#line 151
_stencil_val_a(q2.x,0,0,0);    
    } 
#line 148
{   
      _stencil_val(f,0,0,0);
_stencil_val(u.y,0,0,0);
      
#line 150
_stencil_val_a(q1.y,0,0,0);
_stencil_val(u.y,0,0,0);  
      
#line 151
_stencil_val_a(q2.y,0,0,0);    
    } 
#line 148
{   
      _stencil_val(f,0,0,0);
_stencil_val(u.z,0,0,0);
      
#line 150
_stencil_val_a(q1.z,0,0,0);
_stencil_val(u.z,0,0,0);  
      
#line 151
_stencil_val_a(q2.z,0,0,0);    
    }}end_foreach_stencil()
#line 147 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h"
   BEGIN_FOREACH{foreach()
    { {
      double fc = clamp(val(f,0,0,0),0,1);
      val(q1.x,0,0,0) = fc*rho1*val(u.x,0,0,0);
      val(q2.x,0,0,0) = (1. - fc)*rho2*val(u.x,0,0,0);
    } 
#line 148
{
      double fc = clamp(val(f,0,0,0),0,1);
      val(q1.y,0,0,0) = fc*rho1*val(u.y,0,0,0);
      val(q2.y,0,0,0) = (1. - fc)*rho2*val(u.y,0,0,0);
    } 
#line 148
{
      double fc = clamp(val(f,0,0,0),0,1);
      val(q1.z,0,0,0) = fc*rho1*val(u.z,0,0,0);
      val(q2.z,0,0,0) = (1. - fc)*rho2*val(u.z,0,0,0);
    }}end_foreach();}END_FOREACH 






   {
    _attribute[q2.x.i].inverse = true;
    _attribute[q1.x.i].gradient = _attribute[q2.x.i].gradient = _attribute[u.x.i].gradient;
  } 
#line 159
{
    _attribute[q2.y.i].inverse = true;
    _attribute[q1.y.i].gradient = _attribute[q2.y.i].gradient = _attribute[u.y.i].gradient;
  } 
#line 159
{
    _attribute[q2.z.i].inverse = true;
    _attribute[q1.z.i].gradient = _attribute[q2.z.i].gradient = _attribute[u.z.i].gradient;
  }





  scalar * tracers = _attribute[f.i].tracers;
  _attribute[f.i].tracers = list_concat (tracers, (scalar *)((vector[]){q1, q2,{{-1},{-1},{-1}}}));
  vof_advection (((scalar[]){f,{-1}}), i);
  pfree (_attribute[f.i].tracers,__func__,__FILE__,__LINE__);
  _attribute[f.i].tracers = tracers;





  foreach_stencil(1,{0},)
    {
      {_stencil_val(q1.x,0,0,0); _stencil_val(q2.x,0,0,0);_stencil_val(f,0,0,0);_stencil_val_a(u.x,0,0,0);       }
      
#line 180
{_stencil_val(q1.y,0,0,0); _stencil_val(q2.y,0,0,0);_stencil_val(f,0,0,0);_stencil_val_a(u.y,0,0,0);       }
      
#line 180
{_stencil_val(q1.z,0,0,0); _stencil_val(q2.z,0,0,0);_stencil_val(f,0,0,0);_stencil_val_a(u.z,0,0,0);       }}end_foreach_stencil()





   BEGIN_FOREACH{
#line 178
foreach()
    {
      val(u.x,0,0,0) = (val(q1.x,0,0,0) + val(q2.x,0,0,0))/(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);
      
#line 180
val(u.y,0,0,0) = (val(q1.y,0,0,0) + val(q2.y,0,0,0))/(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);
      
#line 180
val(u.z,0,0,0) = (val(q1.z,0,0,0) + val(q2.z,0,0,0))/(clamp(val(f,0,0,0),0.,1.)*(rho1 - rho2) + rho2);}end_foreach();}END_FOREACH 





  interfaces1 = interfaces, interfaces = NULL;delete((scalar*)((vector[]){q2,q1,{{-1},{-1},{-1}}}));
}{end_tracing("vof_1","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",187);return 0;}end_tracing("vof_1","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",187);}




static int tracer_advection_1_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}





#line 192
      static int tracer_advection_1(const int i,const double t,Event *_ev){tracing("tracer_advection_1","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",192); {
  interfaces = interfaces1;
}{end_tracing("tracer_advection_1","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",194);return 0;}end_tracing("tracer_advection_1","/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",194);}
#line 18 "main3D.c"

#line 1 "view.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/view.h"
#line 56 "/home/esteban/ProgramFile/basilisk/src/view.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/gl/framebuffer.h"
typedef struct _framebuffer framebuffer;

framebuffer * framebuffer_new (unsigned width, unsigned height);
void framebuffer_destroy (framebuffer * p);
unsigned char * framebuffer_image (framebuffer * p);
float * framebuffer_depth (framebuffer * p);
#line 57 "/home/esteban/ProgramFile/basilisk/src/view.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/gl/trackball.h"
#line 50 "/home/esteban/ProgramFile/basilisk/src/gl/trackball.h"
void
gl_trackball(float q[4], float p1x, float p1y, float p2x, float p2y);
#line 61 "/home/esteban/ProgramFile/basilisk/src/gl/trackball.h"
void
gl_add_quats(float q1[4], float q2[4], float dest[4]);





void
gl_build_rotmatrix(float m[4][4], float q[4]);






void
gl_axis_to_quat(float a[3], float phi, float q[4]);
#line 58 "/home/esteban/ProgramFile/basilisk/src/view.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/gl/utils.h"
#line 1 "./gl/tinygl.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/gl/tinygl.h"
#line 122 "/home/esteban/ProgramFile/basilisk/src/gl/tinygl.h"
typedef unsigned int GLenum;
typedef int GLint;
typedef unsigned int GLuint;
typedef float GLfloat;
typedef double GLdouble;
typedef int GLsizei;
typedef unsigned int GLbitfield;
typedef unsigned char GLubyte;

void glBegin (GLenum mode);
void glEnd (void);
void glClear (GLbitfield mask);
void glClearColor (float red, float green, float blue, float alpha);
void glBindTexture (GLenum target, GLuint texture);
void glColor3f (GLfloat red, GLfloat green, GLfloat blue);
void glColorMaterial (GLenum face, GLenum mode);
void glDisable (GLenum cap);
void glEnable (GLenum cap);
void glFinish (void);
void glGetDoublev (GLenum pname, GLdouble * params);
void glGetIntegerv (GLenum pname, GLint * params);
void glHint (GLenum target, GLenum mode);
void glLightfv (GLenum light, GLenum pname, const GLfloat *params);
void glGetLightfv (GLenum light, GLenum pname, GLfloat * params);
void glLightModeli (GLenum pname, GLint param);
void glLineWidth (GLfloat width);
void glPointSize (GLfloat size);
void glNormal3d (GLdouble nx, GLdouble ny, GLdouble nz);
void glOrtho (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top,
       GLdouble nearVal, GLdouble farVal);
void glShadeModel (GLenum mode);
void glTexCoord1d (GLdouble s);
void glTexCoord2f (GLfloat s, GLfloat t);
void glTexImage1D (GLenum target, GLint level, GLint internalFormat,
     GLsizei width, GLint border, GLenum format, GLenum type,
     const void * data);
void glTexParameteri (GLenum target, GLenum pname, GLint param);
void glVertex3d (GLdouble x, GLdouble y, GLdouble z);
void glVertex3f (GLfloat x, GLfloat y, GLfloat z);

GLenum glGetError (void);
void glGetFloatv (GLenum pname, GLfloat * params);
void glMultMatrixf (const GLfloat * m);
void glLoadIdentity (void);
void glScalef (GLfloat x, GLfloat y, GLfloat z);
void glTranslatef (GLfloat x, GLfloat y, GLfloat z);
void glRotatef (GLfloat angle, GLfloat x, GLfloat y, GLfloat z);
void glMatrixMode (GLenum mode);
void glPopMatrix (void);
void glPushMatrix (void);
void glLoadMatrixd (const GLdouble * m);
#line 2 "/home/esteban/ProgramFile/basilisk/src/gl/utils.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 3 "/home/esteban/ProgramFile/basilisk/src/gl/utils.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/stdbool.h"
#include <stdbool.h>
#line 4 "/home/esteban/ProgramFile/basilisk/src/gl/utils.h"

void gl_write_image (FILE * fp, const GLubyte * buffer,
       unsigned width, unsigned height, unsigned samples);
void gl_write_image_png (FILE * fp, const GLubyte * buffer,
    unsigned width, unsigned height, unsigned samples);
void init_gl();

void matrix_multiply (float * m, const float * n);
void vector_multiply (float * v, const float * m);

typedef struct {
  float m[16], p[16];
  float n[6][3];
  float d[6];
  unsigned width;
} Frustum;

void gl_get_frustum (Frustum * f);
int sphere_in_frustum (double x, double y, double z, double r, Frustum * f);
float sphere_diameter (double x, double y, double z, double r, Frustum * f);
void gl_check_error();

int polygonize (const double val[8], double isolevel, double triangles[5][3][3]);
void gl_perspective (double fovy, double aspect, double zNear, double zFar);
int gl_project (float objx, float objy, float objz,
  const float modelMatrix[16],
  const float projMatrix[16],
  const GLint viewport[4],
  float *winx, float *winy, float *winz);

#line 1 "/home/esteban/ProgramFile/basilisk/src/gl/parser.h"
typedef struct _Node Node;

struct _Node {
  char type;
  union {
    char * id;
    double (* func) (double);
    double value;
  } d;
  int s;
  Node * e[3];
};

Node * parse_node (char * code);
void free_node (Node * n);
void print_node (Node * n, FILE * fp);
void reset_node_type (Node * n, char type);
#line 35 "/home/esteban/ProgramFile/basilisk/src/gl/utils.h"
#line 59 "/home/esteban/ProgramFile/basilisk/src/view.h"



#line 1 "input.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/input.h"
#line 16 "/home/esteban/ProgramFile/basilisk/src/input.h"
void input_pgm (scalar s, FILE * fp,
  double ox, double oy, double width)
{
  char line[81];
  if (!fgets (line, 81, fp)) {
    fprintf (ferr, "input_pgm: could not read magic number\n");
    exit (1);
  }
  if (strcmp (line, "P2\n") && strcmp (line, "P5\n")) {
    fprintf (ferr, "input_pgm: magic number '%s' does not match PGM\n",
      line);
    exit (1);
  }
  int binary = !strcmp (line, "P5\n");
  if (!fgets (line, 81, fp)) {
    fprintf (ferr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  int W, H;
  while (line[0] == '#' && fgets (line, 81, fp));
  if (line[0] == '#' || sscanf (line, "%d %d", &W, &H) != 2) {
    fprintf (ferr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  if (!fgets (line, 81, fp)) {
    fprintf (ferr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  int maxval;
  if (sscanf (line, "%d", &maxval) != 1) {
    fprintf (ferr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  if (maxval < 256) {
    unsigned char * a = ((unsigned char *) pmalloc ((W*H)*sizeof(unsigned char),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 1, W*H, fp);
    else {
      int v;
      while (n < W*H && fscanf (fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != W*H) {
      fprintf (ferr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
    foreach_stencil(1,{0},) {          
      
{
 {_stencil_val_a(s,0,0,0);          }
 
{_stencil_val_a(s,0,0,0);  }}
                     
      
    
#line 69
}end_foreach_stencil()
     BEGIN_FOREACH{
#line 63
foreach() {
      int i = (x - ox)*W/width, j = (y - oy)*W/width;
      if (i >= 0 && i < W && j >= 0 && j < H)
 val(s,0,0,0) = 1. - a[(H - 1 - j)*W + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    }end_foreach();}END_FOREACH 
    pfree (a,__func__,__FILE__,__LINE__);
  }
  else {
    unsigned short * a = ((unsigned short *) pmalloc ((W*H)*sizeof(unsigned short),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 2, W*H, fp);
    else {
      int v;
      while (n < W*H && fscanf (fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != W*H) {
      fprintf (ferr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
    foreach_stencil(1,{0},) {          
      
{
 {_stencil_val_a(s,0,0,0);          }
 
{_stencil_val_a(s,0,0,0);  }}
                     
      
    
#line 92
}end_foreach_stencil()
     BEGIN_FOREACH{
#line 86
foreach() {
      int i = (x - ox)*W/width, j = (y - oy)*W/width;
      if (i >= 0 && i < W && j >= 0 && j < H)
 val(s,0,0,0) = 1. - a[(H - 1 - j)*W + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    }end_foreach();}END_FOREACH 
    pfree (a,__func__,__FILE__,__LINE__);
  }
}

static void next_char (FILE * fp, int target)
{
  int c = fgetc(fp), para = 0;
  while (c != EOF && (c != target || para > 0)) {
    if (c == '{') para++;
    if (c == '}') para--;
    c = fgetc(fp);
  }
  if (c != target) {
    fprintf (ferr, "input_gfs(): error: expecting '%c'\n", target);
    exit (1);
  }
}

static int next_string (FILE * fp, const char * target)
{
  int slen = strlen (target), para = 0;
  char s[slen + 1];
  s[slen] = '\0';
  int len = 0, c = fgetc (fp);
  while (c != EOF && len < slen) {
    if (c == '{') para++;
    if (c == '}') para--;
    s[len++] = c;
    c = fgetc (fp);
  }
  while (c != EOF && para >= 0) {
    if (!strcmp (s, target) && para == 0)
      break;
    if (c == '{') para++;
    if (c == '}') para--;
    for (int i = 0; i < slen - 1; i++)
      s[i] = s[i+1];
    s[slen - 1] = c;
    c = fgetc (fp);
  }
  if (strcmp (s, target))
    c = -1;
  return c;
}
#line 156 "/home/esteban/ProgramFile/basilisk/src/input.h"
     
void input_gfs (FILE * fp,
  scalar * list,
  char * file)
{tracing("input_gfs","/home/esteban/ProgramFile/basilisk/src/input.h",157);
  not_mpi_compatible();

  if (file && !(fp = fopen (file, "r"))) {
    perror (file);
    exit (1);
  }

  bool input_all = (list == all);
  if (!list) list = all;





  next_char (fp, '{');

  char * s = ((char *) pmalloc ((1)*sizeof(char),__func__,__FILE__,__LINE__));
  int len = 0;
  int c = fgetc(fp);
  while (c != EOF && c != '}') {
    s[len++] = c;
    s = (char *) prealloc (s, (len + 1)*sizeof(char),__func__,__FILE__,__LINE__);
    s[len] = '\0';
    c = fgetc(fp);
  }
  if (c != '}') {
    fprintf (ferr, "input_gfs(): error: expecting '}'\n");
    exit (1);
  }

  char * s1 = strstr (s, "variables");
  if (!s1) {
    fprintf (ferr, "input_gfs(): error: expecting 'variables'\n");
    exit (1);
  }

  s1 = strstr (s1, "=");
  if (!s1) {
    fprintf (ferr, "input_gfs(): error: expecting '='\n");
    exit (1);
  }
  s1++;

  while (strchr (" \t", *s1))
    s1++;

  scalar * input = NULL;
  s1 = strtok (s1, ", \t");
  while (s1) {
    char * name = replace (s1, '_', '.', false);
    bool found = false;
    {scalar*_i=(scalar*)( list);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
      if (!is_constant(s) && _attribute[s.i].name && !strcmp (_attribute[s.i].name, name)) {
 input = list_append (input, s);
 found = true; break;
      }}}
    if (!found) {
      if (input_all) {
 scalar s = new_scalar("s");
 pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
 _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
 input = list_append (input, s);
      }
      else
 input = list_append (input, (scalar){INT_MAX});
    }
    pfree (name,__func__,__FILE__,__LINE__);
    s1 = strtok (NULL, ", \t");
  }
  pfree (s,__func__,__FILE__,__LINE__);

  next_char (fp, '{');
  double t1 = 0.;
  if (next_string (fp, "Time") >= 0) {
    next_char (fp, '{');
    next_char (fp, 't');
    next_char (fp, '=');
    if (fscanf (fp, "%lf", &t1) != 1) {
      fprintf (ferr, "input_gfs(): error: expecting 't'\n");
      exit (1);
    }
    next_char (fp, '}');
    next_char (fp, '}');
  }

  if (next_string (fp, "Box") < 0) {
    fprintf (ferr, "input_gfs(): error: expecting 'GfsBox'\n");
    exit (1);
  }

  next_char (fp, '{');
  next_char (fp, '{');
  next_char (fp, '\n');

  scalar * listm =((scalar[]) {cm,fm.x,fm.y,fm.z,{-1}});
  scalar * listr = !is_constant(cm) ? listm : NULL;
  NOT_UNUSED (listr);

   BEGIN_FOREACH{foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof (unsigned), 1, fp) != 1) {
      fprintf (ferr, "input_gfs(): error: expecting 'flags'\n");
      exit (1);
    }
    if (!(flags & (1 << 4)) && is_leaf(cell))
      refine_cell (point, listr, 0, NULL);
    double a;
    if (fread (&a, sizeof (double), 1, fp) != 1 || a != -1) {
      fprintf (ferr, "input_gfs(): error: expecting '-1'\n");
      exit (1);
    }
    {scalar*_i=(scalar*)( input);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){ {
      if (fread (&a, sizeof (double), 1, fp) != 1) {
 fprintf (ferr, "input_gfs(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX) {
 if (_attribute[s.i].v.x.i >= 0) {



   if (_attribute[s.i].v.x.i == s.i) {
     s = _attribute[s.i].v.y;
     val(s,0,0,0) = a;
   }
   else if (_attribute[s.i].v.y.i == s.i) {
     s = _attribute[s.i].v.x;
     val(s,0,0,0) = - a;
   }


   else
     val(s,0,0,0) = a;

 }
 else
   val(s,0,0,0) = a;
      }
    }}}
    if (is_leaf(cell))
      continue;
  }end_foreach_cell();}END_FOREACH 
  {scalar*_i=(scalar*)( listm);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s))
      _attribute[s.i].dirty = true;}}
  {scalar*_i=(scalar*)( input);if(_i)for(scalar s=*_i;(&s)->i>=0;s=*++_i){
    if (!is_constant(s))
      _attribute[s.i].dirty = true;}}

  pfree (input,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (t < t1 && events (false))
    t = tnext;
  events (false);
end_tracing("input_gfs","/home/esteban/ProgramFile/basilisk/src/input.h",318);}
#line 357 "/home/esteban/ProgramFile/basilisk/src/input.h"
void input_grd (scalar s,
  FILE * fp, const char * file,
  double nodatavalue,
  bool linear, bool periodic, bool zero,
  int smooth)
{
  scalar input = s;

  if (file && !(fp = fopen (file, "r"))) {
    perror (file);
    exit (1);
  }


  double DeltaGRD;
  int nx, ny;
  double XG0, YG0, ndv;


  char waste[100];
  if (fscanf (fp, "%s %d", waste, &nx) != 2) {
    fprintf (ferr, "input_grd(): error reading 'nx'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %d", waste, &ny) != 2) {
    fprintf (ferr, "input_grd(): error reading 'ny'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %lf", waste, &XG0) != 2) {
    fprintf (ferr, "input_grd(): error reading 'XG0'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %lf", waste, &YG0) != 2) {
    fprintf (ferr, "input_grd(): error reading 'YG0'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %lf", waste, &DeltaGRD) != 2) {
    fprintf (ferr, "input_grd(): error reading 'DeltaGRD'\n");
    if (file) fclose (fp);
    return;
  }
  if (fscanf (fp, "%s %lf", waste, &ndv) != 2) {
    fprintf (ferr, "input_grd(): error reading 'ndv'\n");
    if (file) fclose (fp);
    return;
  }


  if (!nodatavalue)
    nodatavalue = ndv;


  double * value = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,__LINE__));
  for (int i = ny - 1; i >= 0; i--)
    for (int j = 0 ; j < nx; j++) {
      if (fscanf (fp, "%lf ", &value[j + i*nx]) != 1) {
 fprintf (ferr, "input_grd(): error reading value %d,%d\n", i, j);
 if (file) fclose (fp);
 pfree (value,__func__,__FILE__,__LINE__);
 return;
      }
      if (zero && value[j + i*nx] == ndv)
 value[j + i*nx] = 0.;
    }


  if (smooth > 0) {
    double * smoothed = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,__LINE__));
    for (int s = 0; s < smooth; s++) {
      for (int i = 0; i < ny; i++)
 for (int j = 0 ; j < nx; j++) {
   int n = 0;
   smoothed[j + i*nx] = 0.;
   for (int k = -1; k <= 1; k++)
     for (int l = -1; l <= 1; l++)
       if ((l != 0 || k != 0) &&
    i + k >= 0 && i + k < ny &&
    j + l >= 0 && j + l < nx &&
    value[j + l + (i + k)*nx] != ndv)
  smoothed[j + i*nx] += value[j + l + (i + k)*nx], n++;
   if (n == 0)
     smoothed[j + i*nx] = zero ? 0. : ndv;
   else
     smoothed[j + i*nx] /= n;
 }
      do { double * _tmp_ = value; value = smoothed; smoothed = _tmp_; } while(false);
    }
    pfree (smoothed,__func__,__FILE__,__LINE__);
  }

  bool warning = false;
  foreach_stencil (0,{0},) {                   
     
  
           
           
  
     

    
    
{ {              

 
  
          
                 
           
                  
                   
         
                 
      

      
       





{
 {_stencil_val_a(input,0,0,0);  }
 
{_stencil_val_a(input,0,0,0);  }}
                               
                  
              
      
     
         
      
    
#line 486
} 
{
      _stencil_val_a(input,0,0,0);    
      
    }}
                   
    
  
#line 491
}end_foreach_stencil()
  
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
 BEGIN_FOREACH{
#line 452
foreach () {
    if (periodic || _attribute[input.i].boundary[right] == periodic_bc) {
      if (x > XG0 + nx*DeltaGRD)
 x -= nx*DeltaGRD;
      else if (x < XG0)
 x += nx*DeltaGRD;
    }

    int j = (x - XG0 + DeltaGRD/2.)/DeltaGRD;
    int i = (y - YG0 + DeltaGRD/2.)/DeltaGRD;
    if (i >= 0 && i < ny && j >= 0 && j < nx) {
      double val;

      int j1 = (x - XG0)/DeltaGRD;
      int i1 = (y - YG0)/DeltaGRD;
      if (linear && i1 >= 0 && j1 >= 0 && i1 < ny - 1 && j1 < nx - 1 &&
   value[j1 + i1*nx] != ndv && value[j1 + 1 + i1*nx] != ndv &&
   value[j1 + (i1 + 1)*nx] != ndv && value[j1 + 1 + (i1 + 1)*nx] != ndv) {

 double dx = x - (j1*DeltaGRD + XG0);
 double dy = y - (i1*DeltaGRD + YG0);
 val = (value[j1 + i1*nx] +
        dx*(value[j1 + 1 + i1*nx] - value[j1 + i1*nx])/DeltaGRD +
        dy*(value[j1 + (i1 + 1)*nx] - value[j1 + i1*nx])/DeltaGRD +
        dx*dy*(value[j1 + i1*nx] + value[j1 + 1 + (i1 + 1)*nx] -
        value[j1 + (i1 + 1)*nx] - value[j1 + 1 + i1*nx])
        /sq(DeltaGRD));
      }
      else
 val = value[j + i*nx];
      if (val == ndv)
 val(input,0,0,0) = 1e30;
      else
 val(input,0,0,0) = val;
    }
    else {
      val(input,0,0,0) = 1e30;
      warning = true;
    }
  }end_foreach();}END_FOREACH 
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif

  
#line 492
pfree (value,__func__,__FILE__,__LINE__);

  if (warning)
    fprintf (ferr,
      "input_grd(): Warning: Raster data is not covering all"
      " the simulation area\n");

  if (file)
    fclose (fp);
}
#line 63 "/home/esteban/ProgramFile/basilisk/src/view.h"







typedef struct {
  char * expr;
  scalar s;
} cexpr;

static scalar get_cexpr (cexpr * cache, const char * expr)
{
  cexpr * c = cache;
  while (c->expr) {
    if (!strcmp (c->expr, expr)) {


      cexpr tmp = *c;
      while ((c + 1)->expr)
 *c = *(c + 1), c++;
      *c = tmp;
      return c->s;
    }
    c++;
  }
  return (scalar){-1};
}

static cexpr * add_cexpr (cexpr * cache, int maxlen,
     const char * expr, scalar s)
{
  cexpr * c = cache;
  while (c->expr) c++;
  int len = c - cache;
  if (len < maxlen) {
    cache = prealloc (cache, sizeof(cexpr)*(len + 2),__func__,__FILE__,__LINE__);
    c = &cache[len];
  }
  else {

    c = cache;
    pfree (c->expr,__func__,__FILE__,__LINE__);
    scalar s = c->s;
    delete (((scalar[]){s,{-1}}));

    while ((c + 1)->expr)
      *c = *(c + 1), c++;
  }
  c->expr = pstrdup (expr,__func__,__FILE__,__LINE__);
  c->s = s;
  (c + 1)->expr = NULL;
  return cache;
}

static void free_cexpr (cexpr * cache)
{
  cexpr * c = cache;
  while (c->expr) {
    pfree (c->expr,__func__,__FILE__,__LINE__);
    scalar s = c->s;
    delete (((scalar[]){s,{-1}}));
    c++;
  }
  pfree (cache,__func__,__FILE__,__LINE__);
}






typedef void (* MapFunc) (coord *);

struct _bview {
  float tx, ty, sx, sy, sz;
  float quat[4];
  float fov;
  float tz, near, far;

  bool gfsview;
  bool reversed;

  float bg[3];
  float lc;
  float res;

  unsigned width, height, samples;

  framebuffer * fb;
  Frustum frustum;

  MapFunc map;

  int ni;

  bool active;

  cexpr * cache;
  int maxlen;
};

typedef struct _bview bview;




bview * bview_new()
{
  bview * p = ((bview *) pcalloc (1, sizeof(bview),__func__,__FILE__,__LINE__));

  p->tx = p->ty = 0;
  p->sx = p->sy = p->sz = 1.;
  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;
  p->fov = 24.;
  gl_trackball (p->quat, 0.0, 0.0, 0.0, 0.0);




  p->bg[0] = 0.3; p->bg[1] = 0.4; p->bg[2] = 0.6;

  p->res = 1.;
  p->lc = 0.004;

  p->samples = 4;
  p->width = 600*p->samples, p->height = 600*p->samples;


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  p->fb = framebuffer_new (p->width, p->height);

  init_gl();
  p->active = false;

  enable_fpe (FE_DIVBYZERO|FE_INVALID);

  return p;
}




void bview_destroy (bview * p)
{
  framebuffer_destroy (p->fb);
  if (p->cache)
    free_cexpr (p->cache);
  pfree (p,__func__,__FILE__,__LINE__);
}




static bview * _view = NULL;






static void destroy_view()
{
  if (!(_view)) qassert ("/home/esteban/ProgramFile/basilisk/src/view.h", 228, "_view");
  bview_destroy (_view);
}

bview * get_view() {
  if (!_view) {
    _view = bview_new();
    free_solver_func_add (destroy_view);
  }
  return _view;
}




static void redraw() {
  bview * view = get_view();


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  glMatrixMode (0x1701);
  glLoadIdentity ();

  if (view->far <= view->near) {
    double max = 2.;
    gl_perspective (view->fov, view->width/(float)view->height, 1., 1. + 2.*max);

    glMatrixMode (0x1700);
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, - (1. + max));
  }
  else {
    gl_perspective (view->fov, view->width/(float)view->height,
      view->near, view->far);

    glMatrixMode (0x1700);
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, view->tz);
  }

  GLfloat m[4][4];
  gl_build_rotmatrix (m, view->quat);
  glMultMatrixf (&m[0][0]);

  if (view->gfsview) {
    m[0][0] = 0., m[0][1] = 0., m[0][2] = -1.;
    m[1][0] = 0., m[1][1] = -1., m[1][2] = 0.;
    m[2][0] = 1., m[2][1] = 0., m[2][2] = 0.;
    glMultMatrixf (&m[0][0]);
  }

  glScalef (view->sx/L0, view->sy/L0, view->sz/L0);

  glClearColor (view->bg[0], view->bg[1], view->bg[2], 0.);
  glClear (0x00004000|0x00000100);

  gl_get_frustum (&view->frustum);

  view->active = true;
  view->ni = 0;
}




bview * draw() {
  bview * view = get_view();
  if (!view->active)
    redraw();
  else


    disable_fpe (FE_DIVBYZERO|FE_INVALID);
  glMatrixMode (0x1701);
  glTranslatef (0, 0, - 1e-4);
  return view;
}







typedef void * pointer;
#line 357 "/home/esteban/ProgramFile/basilisk/src/view.h"
typedef struct {
  GLubyte a[4];
  float depth;
} RGBA;

static void compose_image_op (void * pin, void * pout, int * len,
         MPI_Datatype * dptr)
{
  RGBA * rin = pin, * out = pout;
  for (int i = 0; i < *len; i++,rin++,out++)
    if (out->depth > rin->depth)
      *out = *rin;
}

     
static pointer compose_image (bview * view)
{tracing("compose_image","/home/esteban/ProgramFile/basilisk/src/view.h",372);
  unsigned char * image = framebuffer_image (view->fb);
  if (!(image)) qassert ("/home/esteban/ProgramFile/basilisk/src/view.h", 375, "image");
  if (npe() > 1) {
    MPI_Op op;
    MPI_Op_create (compose_image_op, true, &op);
    MPI_Datatype rgba;
    MPI_Type_create_struct (2,
       (int[]){4,1},
       (MPI_Aint[]){0,4},
       (MPI_Datatype[]){MPI_BYTE, MPI_FLOAT},
       &rgba);
    MPI_Type_commit (&rgba);
    float * depth = framebuffer_depth (view->fb);
    int size = view->width*view->height;
    RGBA * buf = pmalloc (size*sizeof(RGBA),__func__,__FILE__,__LINE__);
    unsigned char * ptr = image;
    float * dptr = depth;
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < 4; j++)
 buf[i].a[j] = *ptr++;
      buf[i].depth = *dptr++;
    }
    if (pid() == 0) {
      MPI_Reduce (MPI_IN_PLACE, buf, size, rgba, op, 0, MPI_COMM_WORLD);
      unsigned char * ptr = image;
      for (int i = 0; i < size; i++)
 for (int j = 0; j < 4; j++)
   *ptr++ = buf[i].a[j];
    }
    else
      MPI_Reduce (buf, buf, size, rgba, op, 0, MPI_COMM_WORLD);
    pfree (buf,__func__,__FILE__,__LINE__);
    MPI_Op_free (&op);
    MPI_Type_free (&rgba);
  }
  {end_tracing("compose_image","/home/esteban/ProgramFile/basilisk/src/view.h",409);return image;}
end_tracing("compose_image","/home/esteban/ProgramFile/basilisk/src/view.h",410);}



#line 1 "vertexbuffer.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/vertexbuffer.h"
#line 14 "/home/esteban/ProgramFile/basilisk/src/vertexbuffer.h"
struct {

  Array * position, * normal, * color, * index;
  float modelview[16];
  int type;
  int dim;
  int vertex, nvertex;
  bool visible;


  int line_loop, lines, line_strip ;
  int quads, polygon, fan;
  int state;
} VertexBuffer = {
  .visible = false,
  .modelview = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }
};

static void vertex_buffer_push_index (unsigned int i)
{
  i -= VertexBuffer.vertex;
  array_append (VertexBuffer.index, &i, sizeof(unsigned int));
}

void vertex_buffer_setup()
{
  VertexBuffer.nvertex = 0;
  VertexBuffer.type = -1;
  VertexBuffer.dim = -1;
  VertexBuffer.position = array_new();
  VertexBuffer.normal = array_new();
  VertexBuffer.color = array_new();
  VertexBuffer.index = array_new();
}

void vertex_buffer_free()
{
  array_free (VertexBuffer.position);
  VertexBuffer.position = NULL;
  array_free (VertexBuffer.normal);
  VertexBuffer.normal = NULL;
  array_free (VertexBuffer.color);
  VertexBuffer.color = NULL;
  array_free (VertexBuffer.index);
  VertexBuffer.index = NULL;
}

static void vertex_buffer_glBegin (unsigned int state)
{
  if (VertexBuffer.index) {

    glGetFloatv (0x0BA6, VertexBuffer.modelview);

    bview * view = get_view();

    float q[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
      - view->tx, - view->ty, 3, 1 };
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    gl_build_rotmatrix ((float (*)[4])q, view->quat);
    do { float _tmp_ = q[1]; q[1] = q[4]; q[4] = _tmp_; } while(false);
    do { float _tmp_ = q[2]; q[2] = q[8]; q[8] = _tmp_; } while(false);
    do { float _tmp_ = q[6]; q[6] = q[9]; q[9] = _tmp_; } while(false);
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    VertexBuffer.state = state;
    switch (state) {
    case 0x0002:
      VertexBuffer.line_loop = VertexBuffer.nvertex;
      break;
    case 0x0001:
      VertexBuffer.lines = VertexBuffer.nvertex;
      break;
    case 0x0003:
      VertexBuffer.line_strip = VertexBuffer.nvertex;
      break;
    case 0x0007:
      VertexBuffer.quads = VertexBuffer.nvertex;
      break;
    case 0x0009:
      VertexBuffer.polygon = VertexBuffer.nvertex;
      break;
    case 0x0006:
      VertexBuffer.fan = VertexBuffer.nvertex;
      break;
    default:
      fprintf (ferr, "glBegin (%d) not implemented yet\n", state);
      break;
    }
  }
  else
    glBegin (state);
}

static void vertex_buffer_glEnd()
{
  if (VertexBuffer.index) {
    int type = -1;
    switch (VertexBuffer.state) {

    case 0x0002:
      for (int i = VertexBuffer.line_loop; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      vertex_buffer_push_index (VertexBuffer.nvertex - 1);
      vertex_buffer_push_index (VertexBuffer.line_loop);
      type = 0;
      break;

    case 0x0001:
      for (int i = VertexBuffer.lines; i < VertexBuffer.nvertex; i += 2) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;

    case 0x0003:
      for (int i = VertexBuffer.line_strip; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;

    case 0x0007:
      for (int i = VertexBuffer.quads; i < VertexBuffer.nvertex; i += 4)
 for (int j = 1; j <= 2; j++) {
   vertex_buffer_push_index (i);
   vertex_buffer_push_index (i + j);
   vertex_buffer_push_index (i + j + 1);
 }
      type = 1;
      break;

    case 0x0009:
      for (int j = 1; j <= VertexBuffer.nvertex - VertexBuffer.polygon - 2;
    j++) {
 vertex_buffer_push_index (VertexBuffer.polygon);
 vertex_buffer_push_index (VertexBuffer.polygon + j);
 vertex_buffer_push_index (VertexBuffer.polygon + j + 1);
      }
      type = 1;
      break;

    case 0x0006:
      for (int i = VertexBuffer.fan + 1; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (VertexBuffer.fan);
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 1;
      break;

    default:
      break;
    }
    VertexBuffer.state = 0;
    if (VertexBuffer.type >= 0 && type >= 0) {

      if (!(VertexBuffer.type == type)) qassert ("/home/esteban/ProgramFile/basilisk/src/vertexbuffer.h", 179, "VertexBuffer.type == type");
    }
    else
      VertexBuffer.type = type;
  }
  else
    glEnd();
}

static void vertex_buffer_glColor3f (float r, float g, float b)
{
  if (VertexBuffer.color) {
    struct { float x, y, z; } color = {r, g, b};
    array_append (VertexBuffer.color, &color, 3*sizeof(float));
  }
  else
    glColor3f (r, g, b);
}

static void vertex_buffer_glNormal3d (double nx, double ny, double nz)
{
  if (VertexBuffer.normal) {
    struct { float x, y, z; } normal = {nx, ny, nz};
    array_append (VertexBuffer.normal, &normal, 3*sizeof(float));
  }
  else
    glNormal3d (nx, ny, nz);
}

static void vertex_buffer_glVertex3d (double x, double y, double z)
{
  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 3)
      VertexBuffer.dim = 3;
    float v[4] = {x, y, z, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }
  else
    glVertex3d (x, y, z);
}

static void vertex_buffer_glVertex2d (double x, double y)
{
  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 2)
      VertexBuffer.dim = 2;
    float v[4] = {x, y, 0, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }
  else
    glVertex3d (x, y, 0.);
}
#line 415 "/home/esteban/ProgramFile/basilisk/src/view.h"






#line 1 "draw.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/draw.h"





#line 1 "gl/font.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/gl/font.h"
#line 27 "/home/esteban/ProgramFile/basilisk/src/gl/font.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/ast/std/stdio.h"
#include <stdio.h>
#line 28 "/home/esteban/ProgramFile/basilisk/src/gl/font.h"
#line 1 "gl/og_font.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/gl/og_font.h"




#line 1 "gl/tinygl.h"
#line 6 "/home/esteban/ProgramFile/basilisk/src/gl/og_font.h"

typedef struct tagSOG_StrokeVertex SOG_StrokeVertex;
struct tagSOG_StrokeVertex
{
    GLfloat X, Y;
};

typedef struct tagSOG_StrokeStrip SOG_StrokeStrip;
struct tagSOG_StrokeStrip
{
    int Number;
    const SOG_StrokeVertex *Vertices;
};

typedef struct tagSOG_StrokeChar SOG_StrokeChar;
struct tagSOG_StrokeChar
{
    GLfloat Right;
    int Number;
    const SOG_StrokeStrip* Strips;
};

typedef struct tagSOG_StrokeFont SOG_StrokeFont;
struct tagSOG_StrokeFont
{
    char *Name;
    int Quantity;
    GLfloat Height;
    const SOG_StrokeChar **Characters;
};
#line 29 "/home/esteban/ProgramFile/basilisk/src/gl/font.h"
#line 39 "/home/esteban/ProgramFile/basilisk/src/gl/font.h"
extern SOG_StrokeFont ogStrokeMonoRoman;
#line 48 "/home/esteban/ProgramFile/basilisk/src/gl/font.h"
static SOG_StrokeFont *oghStrokeByID( void *font )
{


    if( font == ((void *)0x0001) )
        return &ogStrokeMonoRoman;

    fprintf (ferr, "stroke font %p not found", font );
    return 0;
}
#line 83 "/home/esteban/ProgramFile/basilisk/src/gl/font.h"
void gl_StrokeCharacter( int character )
{
    void *fontID = ((void *)0x0001);
    const SOG_StrokeChar *schar;
    const SOG_StrokeStrip *strip;
    int i, j;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( !font ||
        ( 1 > character ) ||
        ( font->Quantity < character ) )
        return;

    schar = font->Characters[ character ];
    if( schar )
    {
        strip = schar->Strips;

        for( i = 0; i < schar->Number; i++, strip++ )
        {
            vertex_buffer_glBegin( 0x0003 );
            for( j = 0; j < strip->Number; j++ )
                vertex_buffer_glVertex2d( strip->Vertices[ j ].X, strip->Vertices[ j ].Y );
            vertex_buffer_glEnd( );
        }
        glTranslatef( schar->Right, 0.0, 0.0 );
    }
}
#line 147 "/home/esteban/ProgramFile/basilisk/src/gl/font.h"
void gl_StrokeString( const char *string )
{
    void *fontID = ((void *)0x0001);
    int i, j;
    float length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );
    unsigned char c;

    if( font && string )





        while(( c = *string++ ))
       if( c < font->Quantity ) {
                if( c == '\n' )
                {
                    glTranslatef ( -length, -( float )( font->Height ), 0.0 );
                    length = 0.0;
                }
                else
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                    {
                        const SOG_StrokeStrip *strip = schar->Strips;

                        for( i = 0; i < schar->Number; i++, strip++ )
                        {
                            vertex_buffer_glBegin( 0x0003 );

                            for( j = 0; j < strip->Number; j++ )
                                vertex_buffer_glVertex2d( strip->Vertices[ j ].X,
                                            strip->Vertices[ j ].Y);

                            vertex_buffer_glEnd( );
                        }

                        length += schar->Right;
                        glTranslatef( schar->Right, 0.0, 0.0 );
                    }
                }
     }
}
#line 226 "/home/esteban/ProgramFile/basilisk/src/gl/font.h"
float gl_StrokeWidth( int character )
{
    void *fontID = ((void *)0x0001);
    float ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font &&
        ( 0 < character ) &&
        ( font->Quantity > character ) )
    {
        const SOG_StrokeChar *schar = font->Characters[ character ];
        if( schar )
            ret = schar->Right;
    }

    return ret;
}
#line 269 "/home/esteban/ProgramFile/basilisk/src/gl/font.h"
float gl_StrokeLength( const char *string )
{
    void *fontID = ((void *)0x0001);
    unsigned char c;
    float length = 0.0;
    float this_line_length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font && string )
        while(( c = *string++ ))
            if( c < font->Quantity )
            {
                if( c == '\n' )
                {
                    if( length < this_line_length )
                        length = this_line_length;
                    this_line_length = 0.0;
                }
                else
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                        this_line_length += schar->Right;
                }
            }

    if( length < this_line_length )
        length = this_line_length;
    return length;
}
#line 321 "/home/esteban/ProgramFile/basilisk/src/gl/font.h"
GLfloat gl_StrokeHeight()
{
    void *fontID = ((void *)0x0001);
    GLfloat ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font )
        ret = font->Height;

    return ret;
}
#line 7 "/home/esteban/ProgramFile/basilisk/src/draw.h"




void clear()
{
  bview * view = get_view();
  if (view->active)
    view->active = false;
  draw();
}
#line 49 "/home/esteban/ProgramFile/basilisk/src/draw.h"
void view (float tx, float ty,
    float fov,
    float quat[4],
    float sx, float sy, float sz,
    unsigned width, unsigned height, unsigned samples,
    float bg[3],
    float theta, float phi, float psi,
    bool relative,
    float tz, float near, float far,
    float res,
    char * camera,
    MapFunc map,
    int cache,
    float p1x, float p1y, float p2x, float p2y,
    bview * view1)
{
  bview * v = view1 ? view1 : get_view();
  if (fov) {
    if (relative)
      v->fov += (0.1 + 3.*v->fov)*fov;
    else
      v->fov = fov;
    v->fov = clamp(v->fov,0.01,100.);
  }
  for (int i = 0; i < 4; i++)
    if (quat[i]) {
      for (int j = 0; j < 4; j++)
 v->quat[j] = quat[j];
      break;
    }
  v->tx = relative ? v->tx + tx*0.02*(0.01 + 3.*v->fov) : tx;
  v->ty = relative ? v->ty + ty*0.02*(0.01 + 3.*v->fov) : ty;
  v->sx = sx;
  v->sy = sy;
  v->sz = sz;
  if (bg[0] || bg[1] || bg[2])
    for (int i = 0; i < 3; i++)
      v->bg[i] = bg[i];

  if (camera) {
    v->gfsview = false;
    if (strlen(camera) >= 4 &&
 !strcmp (&camera[strlen(camera) - 4], ".gfv")) {
      FILE * fp = fopen (camera, "r");
      if (!fp) {
 perror (camera);
 exit (1);
      }
      char s[81];
      float q[4], fov;
      int nq = 0, nf = 0;
      while (fgets (s, 81, fp) && (!nq || !nf)) {
 if (!nq)
   nq = sscanf (s, "  q0 = %f q1 = %f q2 = %f q3 = %f",
         &q[0], &q[1], &q[2], &q[3]);
 if (!nf)
   nf = sscanf (s, "  fov = %f", &fov);
      }
      if (nq != 4 || nf != 1) {
 fprintf (ferr, "%s: not a valid gfv file\n", camera);
 exit (1);
      }
      for (int j = 0; j < 4; j++)
 v->quat[j] = q[j];
      v->fov = fov;
      v->gfsview = true;
    }
    else if (!strcmp (camera, "left"))
      gl_axis_to_quat ((float[]){0,1,0}, - 3.14159265358979/2., v->quat);
    else if (!strcmp (camera, "right"))
      gl_axis_to_quat ((float[]){0,1,0}, 3.14159265358979/2., v->quat);
    else if (!strcmp (camera, "top"))
      gl_axis_to_quat ((float[]){1,0,0}, - 3.14159265358979/2., v->quat);
    else if (!strcmp (camera, "bottom"))
      gl_axis_to_quat ((float[]){1,0,0}, 3.14159265358979/2., v->quat);
    else if (!strcmp (camera, "front"))
      gl_axis_to_quat ((float[]){0,0,1}, 0., v->quat);
    else if (!strcmp (camera, "back"))
      gl_axis_to_quat ((float[]){0,1,0}, 3.14159265358979, v->quat);
    else if (!strcmp (camera, "iso")) {
      gl_axis_to_quat ((float[]){0,1,0}, 3.14159265358979/4., v->quat);
      float q[4];
      gl_axis_to_quat ((float[]){1,0,0}, - 3.14159265358979/4., q);
      gl_add_quats(q, v->quat, v->quat);
    }
    else {
      fprintf (ferr, "view(): unknown camera '%s'\n", camera);
      exit (1);
    }
  }
  else if (theta || phi || psi) {
    v->gfsview = false;
    float q[4];
    gl_axis_to_quat ((float[]){1,0,0}, - phi, q);
    if (relative) {
      float q1[4];
      gl_axis_to_quat ((float[]){0,1,0}, theta, q1);
      gl_add_quats(q, q1, q1);
      float q2[4];
      gl_axis_to_quat ((float[]){0,0,1}, psi, q2);
      gl_add_quats(q1, q2, q2);
      gl_add_quats(q2, v->quat, v->quat);
    }
    else {
      gl_axis_to_quat ((float[]){0,1,0}, theta, v->quat);
      gl_add_quats(q, v->quat, v->quat);
      gl_axis_to_quat ((float[]){0,0,1}, psi, q);
      gl_add_quats(q, v->quat, v->quat);
    }
  }

  if (map)
    v->map = map;

  if (p1x || p1y || p2x || p2y) {
    float q[4];
    gl_trackball(q, p1x, p1y, p2x, p2y);
    gl_add_quats (q, v->quat, v->quat);
  }

  if (far > near) {
    v->tz = tz;
    v->far = far;
    v->near = near;
  }

  if (res)
    v->res = res;

  if ((width && width != v->width) ||
      (height && height != v->height) ||
      (samples && samples != v->samples)) {
    v->width = v->width/v->samples;
    v->height = v->height/v->samples;
    if (width) v->width = width;
    if (height) v->height = height;
    if (samples) v->samples = samples;
    v->width *= v->samples;
    v->height *= v->samples;
    framebuffer_destroy (v->fb);


    disable_fpe (FE_DIVBYZERO|FE_INVALID);

    v->fb = framebuffer_new (v->width, v->height);
    init_gl();

    enable_fpe (FE_DIVBYZERO|FE_INVALID);
  }

  if (cache > 0) {
    v->cache = pcalloc (1, sizeof (cexpr),__func__,__FILE__,__LINE__);
    v->maxlen = cache;
  }

  clear();
}







void begin_translate (float x, float y, float z)
{
  bview * view = draw();
  glMatrixMode (0x1700);
  glPushMatrix();
  glTranslatef (x, y, z);
  gl_get_frustum (&view->frustum);
}

void end_translate()
{
  bview * view = draw();
  glMatrixMode (0x1700);
  glPopMatrix();
  gl_get_frustum (&view->frustum);
}
#line 238 "/home/esteban/ProgramFile/basilisk/src/draw.h"
void begin_mirror (coord n, double alpha)
{
  bview * view = draw();
  glMatrixMode (0x1700);
  glPushMatrix();
  normalize (&n);
  GLfloat s[16], t[16];
  s[0] = 1. - 2.*n.x*n.x;
  s[1] = - 2.*n.x*n.y; s[2] = - 2.*n.x*n.z;
  s[3] = 0.;
  s[4] = s[1];
  s[5] = 1. - 2.*n.y*n.y; s[6] = - 2.*n.y*n.z;
  s[7] = 0.;
  s[8] = s[2]; s[9] = s[6]; s[10] = 1. - 2.*n.z*n.z;
  s[11] = 0.;
  s[12] = 0.; s[13] = 0.; s[14] = 0.;
  s[15] = 1.;

  t[0] = 1.; t[1] = 0.; t[2] = 0.; t[3] = 0.;
  t[4] = 0.; t[5] = 1.; t[6] = 0.; t[7] = 0.;
  t[8] = 0.; t[9] = 0.; t[10] = 1.; t[11] = 0.;
  t[12] = - 2.*n.x*alpha;
  t[13] = - 2.*n.y*alpha;
  t[14] = - 2.*n.z*alpha;
  t[15] = 1.;
  matrix_multiply (s, t);
  glMultMatrixf (s);
  gl_get_frustum (&view->frustum);
  view->reversed = !view->reversed;
}

void end_mirror() {
  end_translate();
  bview * view = draw();
  view->reversed = !view->reversed;
}







static void mapped_position (bview * view, coord * p, double * r)
{
  double x = p->x, y = p->y, z = p->z, rm = 0.;
  view->map (p);
  for (int i = -1; i <= 1; i += 2)
    for (int j = -1; j <= 1; j += 2)
      for (int k = -1; k <= 1; k += 2) {
 coord q = {x + i**r, y + j**r, z + k**r};
 view->map (&q);
 double pq = sq(p->x - q.x) + sq(p->y - q.y) + sq(p->z - q.z);
 if (pq > rm)
   rm = pq;
      }
  *r = sqrt (rm);
}

#define foreach_visible(view)\
foreach_cell() {\
\
\
\
  double _r = Delta*0.87;\
\
  coord _p = {x, y, z};\
  if ((view)->map)\
    mapped_position (view, &_p, &_r);\
  if (VertexBuffer.visible &&\
      !sphere_in_frustum (_p.x, _p.y, _p.z, _r, &(view)->frustum))\
    continue;\
  if (is_leaf(cell) ||\
      (VertexBuffer.visible &&\
       sphere_diameter (_p.x, _p.y, _p.z, _r/L0, &(view)->frustum)\
       < (view)->res)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 315

#define end_foreach_visible()\
    }\
    continue;\
  }\
}\
end_foreach_cell();\

#line 322

#line 333 "/home/esteban/ProgramFile/basilisk/src/draw.h"
static void glnormal3d (bview * view, double x, double y, double z) {

  if (view->gfsview || view->reversed)
    vertex_buffer_glNormal3d (- x, - y, - z);
  else
    vertex_buffer_glNormal3d (x, y, z);
}

#define foreach_visible_plane(view, n1, alpha1)\
coord _n = {(n1).x, (n1).y, (n1).z};\
double _alpha = 0.9999999*(alpha1);\
{\
  double norm = sqrt(sq(_n.x) + sq(_n.y) + sq(_n.z));\
  if (!norm)\
    _n.z = 1.;\
  else\
    _n.x /= norm, _n.y /= norm, _n.z /= norm, _alpha /= norm;\
}\
glnormal3d (view, _n.x, _n.y, _n.z);\
foreach_cell() {\
\
  double _r = Delta*0.87, alpha = (_alpha - _n.x*x - _n.y*y - _n.z*z)/Delta;\
  if (fabs(alpha) > 0.87 ||\
      (VertexBuffer.visible &&\
       !sphere_in_frustum (x, y, z, _r, &(view)->frustum)))\
    continue;\
  if (is_leaf(cell) ||\
      (VertexBuffer.visible &&\
       sphere_diameter (x, y, z, _r/L0, &(view)->frustum) < (view)->res)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 363

#define end_foreach_visible_plane()\
    }\
    continue;\
  }\
}\
end_foreach_cell();\

#line 370



static bool _reversed = false;

static void begin_draw_lines (bview * view, float color[3], float lw)
{
  glMatrixMode (0x1701);
  glPushMatrix();
  glTranslatef (0., 0., view->lc*view->fov/24.);
  vertex_buffer_glColor3f (color[0], color[1], color[2]);
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));
  _reversed = view->reversed;
  view->reversed = false;
}

static void end_draw_lines()
{
  glMatrixMode (0x1701);
  glPopMatrix();
  bview * view = draw();
  view->reversed = _reversed;
}

static inline double interp (Point point, coord p, scalar col) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  return interpolate_linear (point, col
,
        
#line 396
x + p.x*Delta, y + p.y*Delta, z + p.z*Delta);
}


#line 394
static void _stencil_interp (Point point,_stencil_undefined * p, scalar col) {int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES; 
_stencil_interpolate_linear (point, col
,NULL,NULL,NULL   );
  
#line 395
return
;
}

static double evaluate_expression (Point point, Node * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  if (!(n)) qassert ("/home/esteban/ProgramFile/basilisk/src/draw.h", 401, "n");
  switch (n->type) {
  case '1': return n->d.value;
  case '+': return (evaluate_expression (point, n->e[0]) +
      evaluate_expression(point, n->e[1]));
  case '-': return (evaluate_expression (point, n->e[0]) -
      evaluate_expression(point, n->e[1]));
  case '*': return (evaluate_expression (point, n->e[0]) *
      evaluate_expression(point, n->e[1]));
  case '/': return (evaluate_expression (point, n->e[0]) /
      evaluate_expression(point, n->e[1]));
  case '^': return pow (evaluate_expression (point, n->e[0]),
   evaluate_expression(point, n->e[1]));
  case '>': return (evaluate_expression (point, n->e[0]) >
      evaluate_expression(point, n->e[1]));
  case '<': return (evaluate_expression (point, n->e[0]) <
      evaluate_expression(point, n->e[1]));
  case 'L': return (evaluate_expression (point, n->e[0]) <=
      evaluate_expression(point, n->e[1]));
  case 'G': return (evaluate_expression (point, n->e[0]) >=
      evaluate_expression(point, n->e[1]));
  case '=': return (evaluate_expression (point, n->e[0]) ==
      evaluate_expression(point, n->e[1]));
  case 'i': return (evaluate_expression (point, n->e[0]) !=
      evaluate_expression(point, n->e[1]));
  case 'O': return (evaluate_expression (point, n->e[0]) ||
      evaluate_expression(point, n->e[1]));
  case 'A': return (evaluate_expression (point, n->e[0]) &&
      evaluate_expression(point, n->e[1]));
  case '?': return (evaluate_expression (point, n->e[0]) ?
      evaluate_expression(point, n->e[1]) :
      evaluate_expression(point, n->e[2]));
  case 'm': return - evaluate_expression (point, n->e[0]);
  case 'f': return n->d.func (evaluate_expression (point, n->e[0]));
  case 'v': {
    scalar s = {n->s};
    int k[3] = {0,0,0};
    for (int i = 0; i < 3; i++)
      if (n->e[i])
 k[i] = evaluate_expression (point, n->e[i]);
    return val(s,k[0],k[1],k[2]);
  }
  case 'D': return Delta;
  case 'x': return x;
  case 'y': return y;
  case 'z': return z;
  default:
    fprintf (ferr, "unknown operation type '%c'\n", n->type);
    if (!(false)) qassert ("/home/esteban/ProgramFile/basilisk/src/draw.h", 449, "false");
  }
  return undefined;
}


#line 399
static void _stencil_evaluate_expression (Point point, Node * n)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES; 
      
  switch (n->type) {
  case '1': return ;
  case '+': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 404
return  
;}
  case '-': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 406
return  
;}
  case '*': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 408
return  
;}
  case '/': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 410
return  
;}
  case '^': {_stencil_evaluate_expression (point, n->e[0]);
   _stencil_evaluate_expression(point, n->e[1]);
#line 412
return  
;}
  case '>': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 414
return  
;}
  case '<': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 416
return  
;}
  case 'L': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 418
return  
;}
  case 'G': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 420
return  
;}
  case '=': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 422
return  
;}
  case 'i': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 424
return  
;}
  case 'O': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 426
return  
;}
  case 'A': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
#line 428
return  
;}
  case '?': {_stencil_evaluate_expression (point, n->e[0]);
      _stencil_evaluate_expression(point, n->e[1]);
      _stencil_evaluate_expression(point, n->e[2]);
#line 430
return   

;}
  case 'm': { _stencil_evaluate_expression (point, n->e[0]);return ;}
  case 'f': {_stencil_evaluate_expression (point, n->e[0]);return  ;}
  case 'v': {
    scalar s = {n->s};   
    
    for (int i = 0; i < 3; i++)
      if (n->e[i])
 { _stencil_evaluate_expression (point, n->e[i]); } 
_stencil_val(s,o_stencil,o_stencil,o_stencil);
    
#line 441
return;
  }
  case 'D': return ;
  case 'x': return ;
  case 'y': return ;
  case 'z': return ; 
     
    
        
  }
  return ;
}

static bool assemble_node (Node * n)
{
  if (n->type == 'v') {
    char * id = n->d.id;
    scalar s = lookup_field (id);
    if (s.i >= 0)
      n->s = s.i;
    else {
      n->s = -1;
      if (!strcmp (id, "Delta"))
 reset_node_type (n, 'D');
      else if (!strcmp (id, "x"))
 reset_node_type (n, 'x');
      else if (!strcmp (id, "y"))
 reset_node_type (n, 'y');
      else if (!strcmp (id, "z"))
 reset_node_type (n, 'z');
      else {
 typedef struct { char * name; double val; } Constant;
 static Constant constants[] = {
   {"pi", 3.14159265358979 },
   {"nodata", 1e30 },
   {"HUGE", 1e30 },
   { NULL },
 };
 Constant * p = constants;
 while (p->name) {
   if (!strcmp (p->name, id)) {
     reset_node_type (n, '1');
     n->d.value = p->val;
     break;
   }
   p++;
 }
 if (n->type == 'v') {
   fprintf (ferr, "unknown identifier '%s'\n", id);
   return false;
 }
      }
    }
  }
  for (int i = 0; i < 3; i++)
    if (n->e[i] && !assemble_node (n->e[i]))
      return false;
  return true;
}

static scalar compile_expression (char * expr, bool * isexpr)
{
  *isexpr = false;
  if (!expr)
    return (scalar){-1};

  bview * view = get_view();
  scalar s;
  if (view->cache && (s = get_cexpr (view->cache, expr)).i >= 0)
    return s;

  Node * node = parse_node (expr);
  if (node == NULL) {
    fprintf (ferr, "'%s': syntax error\n", expr);
    return (scalar){-1};
  }
  if (!assemble_node (node)) {
    free_node (node);
    return (scalar){-1};
  }
  if (node->type == 'v' && node->e[0] == NULL) {
    scalar s = {node->s};
    if (_attribute[s.i].block > 0) {
      free_node (node);
      return s;
    }
  }
  s = new_scalar("s");
  pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
  _attribute[s.i].name = pstrdup (expr,__func__,__FILE__,__LINE__);
  foreach_stencil(1,{0},)
    { _stencil_evaluate_expression (point, node);_stencil_val_a(s,0,0,0); }end_foreach_stencil()
   BEGIN_FOREACH{
#line 531
foreach()
    val(s,0,0,0) = evaluate_expression (point, node);end_foreach();}END_FOREACH 
  restriction (((scalar[]){s,{-1}}));
  free_node (node);

  if (view->cache)
    view->cache = add_cexpr (view->cache, view->maxlen, expr, s);
  else
    *isexpr = true;
  return s;
}
#line 604 "/home/esteban/ProgramFile/basilisk/src/draw.h"
static void begin_colorized (float fc[3], bool constant_color,
        double cmap[127][3], bool use_texture)
{

  if (use_texture) {
    GLfloat texture[3*256];
    for (int i = 0; i < 256; i++) {
      Color j = colormap_color (cmap, i/255., 0, 1);
      texture[3*i] = j.r/255.;
      texture[3*i + 1] = j.g/255.;
      texture[3*i + 2] = j.b/255.;
    }
    glTexImage1D (0x0DE0, 0, 0x1907, 256,0, 0x1907, 0x1406, texture);
    glTexParameteri (0x0DE0, 0x2801, 0x2601);
    glTexParameteri (0x0DE0, 0x2800, 0x2601);
    glTexParameteri (0x0DE0, 0x2802, 0x812F);
    glTexParameteri (0x0DE0, 0x2803, 0x812F);
    glEnable (0x0DE0);
  }
  if (constant_color)
    vertex_buffer_glColor3f (fc[0], fc[1], fc[2]);
}

static void end_colorized() {
  glDisable (0x0DE0);
}
#line 656 "/home/esteban/ProgramFile/basilisk/src/draw.h"
     
bool colorbar (Colormap map, float size, float pos[2],
        char * label, double lscale, double min,
        double max, bool horizontal, bool border,
        bool mid, float lc[3], float lw, float fsize,
        char * format, int levels)
{tracing("colorbar","/home/esteban/ProgramFile/basilisk/src/draw.h",657);
  bview * view = draw();
  glDisable (0x0B50);
  glMatrixMode (0x1701);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode (0x1700);
  glPushMatrix();
  glLoadIdentity();

  float fheight = gl_StrokeHeight();
  if (!size)
    size = 15;
  float width = 2./size;
  if (levels < 1) levels = 1;
  float h = 0, height = 4*width, dh = height/levels;
  glTranslatef (pos[0], pos[1], 0);


  double cmap [127][3];
  (* map) (cmap);
  vertex_buffer_glBegin(0x0007);
  for (int i = 0; i < levels; i++) {
    Color c = colormap_color (cmap, (float)i/(levels - 1), 0, 1);
    vertex_buffer_glColor3f (c.r/255., c.g/255., c.b/255.);
    if (horizontal) {
      vertex_buffer_glVertex2d (h + dh, 0);
      vertex_buffer_glVertex2d (h + dh, width);
      vertex_buffer_glVertex2d (h, width);
      vertex_buffer_glVertex2d (h, 0);
    } else {
      vertex_buffer_glVertex2d (0, h + dh);
      vertex_buffer_glVertex2d (width, h + dh);
      vertex_buffer_glVertex2d (width, h);
      vertex_buffer_glVertex2d (0, h);
    }
    h += dh;
    view->ni++;
  }
  vertex_buffer_glEnd();
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));
  vertex_buffer_glColor3f (lc[0], lc[1], lc[2]);


  if (border) {
    vertex_buffer_glBegin (0x0002);
    vertex_buffer_glVertex2d (0, 0);
    if (horizontal) {
      vertex_buffer_glVertex2d (0, width);
      vertex_buffer_glVertex2d (height, width);
      vertex_buffer_glVertex2d (height, 0);
    } else {
      vertex_buffer_glVertex2d (width, 0);
      vertex_buffer_glVertex2d (width, height);
      vertex_buffer_glVertex2d (0, height);
    }
    vertex_buffer_glEnd();
  }


  float fwidth = gl_StrokeWidth ('1');
  if (!fsize)
    fsize = 20;
  float hscale = 2./(fsize*fwidth), vscale = hscale*view->width/view->height;
  char str[99];
  vertex_buffer_glColor3f (lc[0], lc[1], lc[2]);
  if (horizontal)
    glTranslatef (0, -(fheight/(view->height)), 0);
  else
    glTranslatef (width, -(fheight/(3*view->height)), 0);
  glScalef (hscale, vscale, 1.);
  sprintf (str, format, min);
  if (min > -1e30) {
    glPushMatrix();
    if (horizontal)
      glTranslatef (-fwidth*(strlen(str) - 1)/2, 0, 0);
    glScalef (lscale, lscale, 1.);
    gl_StrokeString (str);
    glPopMatrix();
  }
  if (horizontal)
    glTranslatef (height/hscale,0, 0);
  else
    glTranslatef (0, height/vscale, 0);
  sprintf (str, format, max);
  if (max < 1e30) {
    glPushMatrix();
    if (horizontal)
      glTranslatef (-fwidth*(strlen(str) - 1)/2, 0, 0);
    glScalef (lscale, lscale, 1.);
    gl_StrokeString (str);
    glPopMatrix();
  }

  if (mid) {
    sprintf (str, format, (min + max)/2);
    glPushMatrix();
    if (horizontal)
      glTranslatef (-height/(2*hscale) - fwidth*(strlen(str) - 1)/2,0, 0);
    else
      glTranslatef (0, -height/(2*vscale), 0);
    glScalef (lscale, lscale, 1.);
    gl_StrokeString (str);
    glPopMatrix();
  }

  if (horizontal)
    glTranslatef (-height/(2*hscale) - lscale*fwidth*(strlen(label) - 1)/2, width/vscale, 0);
  else
    glTranslatef (-width/hscale, 0, 0);

  glScalef (lscale, lscale, 1.);
  glTranslatef (0, fheight, 0);
  gl_StrokeString (label);

  glMatrixMode (0x1700);
  glPopMatrix();
  glMatrixMode (0x1701);
  glPopMatrix();
  {end_tracing("colorbar","/home/esteban/ProgramFile/basilisk/src/draw.h",781);return true;}
end_tracing("colorbar","/home/esteban/ProgramFile/basilisk/src/draw.h",782);}
#line 801 "/home/esteban/ProgramFile/basilisk/src/draw.h"
static bool cfilter (Point point, scalar c, double cmin)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;
  double cmin1 = 4.*cmin;
  if (val(c,0,0,0) <= cmin) {
    
      if (val(c,1,0,0) >= 1. - cmin1 || val(c,-1,0,0) >= 1. - cmin1)
 return true;
      
#line 806
if (val(c,0,1,0) >= 1. - cmin1 || val(c,0,-1,0) >= 1. - cmin1)
 return true;
      
#line 806
if (val(c,0,0,1) >= 1. - cmin1 || val(c,0,0,-1) >= 1. - cmin1)
 return true;
    return false;
  }
  if (val(c,0,0,0) >= 1. - cmin) {
    
      if (val(c,1,0,0) <= cmin1 || val(c,-1,0,0) <= cmin1)
 return true;
      
#line 812
if (val(c,0,1,0) <= cmin1 || val(c,0,-1,0) <= cmin1)
 return true;
      
#line 812
if (val(c,0,0,1) <= cmin1 || val(c,0,0,-1) <= cmin1)
 return true;
    return false;
  }
  int n = 0;
  double min = 1e30, max = - 1e30;
  {foreach_neighbor(1) {
    if (val(c,0,0,0) > cmin && val(c,0,0,0) < 1. - cmin && ++n >= (1 << 3))
      return true;
    if (val(c,0,0,0) > max) max = val(c,0,0,0);
    if (val(c,0,0,0) < min) min = val(c,0,0,0);
  }end_foreach_neighbor()}
  return max - min > 0.5;
}
#line 801 "/home/esteban/ProgramFile/basilisk/src/draw.h"
static void _stencil_cfilter (Point point, scalar c,_stencil_undefined * cmin)
{int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;   
  
_stencil_val(c,0,0,0); {
    
      {_stencil_val(c,1,0,0); _stencil_val(c,-1,0,0); 
           }
      
#line 806
{_stencil_val(c,0,1,0); _stencil_val(c,0,-1,0); 
           }
      
#line 806
{_stencil_val(c,0,0,1); _stencil_val(c,0,0,-1); 
           } 
    
  }
_stencil_val(c,0,0,0); {
    
      {_stencil_val(c,1,0,0); _stencil_val(c,-1,0,0); 
       }
      
#line 812
{_stencil_val(c,0,1,0); _stencil_val(c,0,-1,0); 
       }
      
#line 812
{_stencil_val(c,0,0,1); _stencil_val(c,0,0,-1); 
       } 
    
  }          
     
       
  
  
  {
#line 818
foreach_neighbor(1) {
_stencil_val(c,0,0,0); _stencil_val(c,0,0,0);

_stencil_val(c,0,0,0); { _stencil_val(c,0,0,0); }
_stencil_val(c,0,0,0); { _stencil_val(c,0,0,0); } 
      
                  
       
       
  
#line 823
}end_foreach_neighbor()}
  return     ;
}

static void glvertex3d (bview * view, double x, double y, double z) {
  if (view->map) {
    coord p = {x, y, z};
    view->map (&p);
    vertex_buffer_glVertex3d (p.x, p.y, p.z);
  }
  else
    vertex_buffer_glVertex3d (x, y, z);
}
#line 885 "/home/esteban/ProgramFile/basilisk/src/draw.h"
     
bool draw_vof (char * c, char * s, bool edges,
        double larger, int filled,
        char * color,
        double min, double max, double spread,
        bool linear,
        Colormap map,
        float fc[3], float lc[3], float lw,
        bool expr,
        bool cbar, float size, float pos[2], char * label, double lscale, bool horizontal, bool border, bool mid, float fsize, char * format, int levels)
{tracing("draw_vof","/home/esteban/ProgramFile/basilisk/src/draw.h",886);
  scalar d = lookup_field (c);
  if (d.i < 0) {
    fprintf (ferr, "draw_vof(): no field named '%s'\n", c);
    {end_tracing("draw_vof","/home/esteban/ProgramFile/basilisk/src/draw.h",899);return false;}
  }
  vector fs = lookup_vector (s);

  scalar col = {-1}; if (color && strcmp (color, "level")) { col = compile_expression (color, &expr); if (col.i < 0) {end_tracing("draw_vof","/home/esteban/ProgramFile/basilisk/src/draw.h",903);return false;} } double cmap[127][3]; if (color) { if (min == 0 && max == 0) { if (col.i < 0) min = 0, max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (spread < 0.) min = s.min, max = s.max; else { if (!spread) spread = 5.; min = avg - spread*s.stddev; max = avg + spread*s.stddev; } } } if (!map) map = jet; (* map) (cmap); } if ((3 > 2 || linear) && !fc[0] && !fc[1] && !fc[2]) fc[0] = fc[1] = fc[2] = 1.; if (cbar) colorbar (map, size, pos, label, lscale, min, max, horizontal, border, mid, lc, lw, fsize, format, levels);;

  double cmin = 1e-3;
#line 916 "/home/esteban/ProgramFile/basilisk/src/draw.h"
  bview * view = draw();
#line 994 "/home/esteban/ProgramFile/basilisk/src/draw.h"
  if (!larger)
    larger = edges || (color && !linear) ? 1. : 1.1;
  if (edges)
    {begin_draw_lines (view, lc, lw); {
      foreach_visible_stencil (1,DEPAREN({view}),)
 {_stencil_cfilter (point, d, NULL); {  
    _stencil_facet_normal (point, d, fs);     
   _stencil_val(d,0,0,0);        
   
      
             
     
             
      
       
     
     
    
 } }end_foreach_visible_stencil()
       BEGIN_FOREACH{
#line 998
foreach_visible (view)
 if (cfilter (point, d, cmin)) {
   coord n = facet_normal (point, d, fs);
   double alpha = plane_alpha (val(d,0,0,0), n);
   coord v[12];
   int m = facets (n, alpha, v, larger);
   if (m > 2) {
     vertex_buffer_glBegin (0x0002);
     for (int i = 0; i < m; i++)
       glvertex3d (view,
     x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
     vertex_buffer_glEnd ();
     view->ni++;
   }
 }end_foreach_visible();}END_FOREACH 
    }end_draw_lines();}
  else
    {begin_colorized (fc, !VertexBuffer.color || !color, cmap, !VertexBuffer.color && color && linear && col.i >= 0); {
      foreach_visible_stencil (1,DEPAREN({view}),)
 {_stencil_cfilter (point, d, NULL); {  
    _stencil_facet_normal (point, d, fs);     
   _stencil_val(d,0,0,0);        
   
    
{    
     
           {
       if (linear) {
  if (color && linear && col.i >= 0) { if (VertexBuffer.color) {        _stencil_interp (point,NULL , col);     } else {    _stencil_interp (point,NULL , col);                } }
       }
       else {
  if (color && (!linear || col.i < 0)) {               _stencil_val(col,0,0,0);     }
       }        
          
       
       
     } 
     
     
   }
      
 
#line 1038
} }end_foreach_visible_stencil()
       BEGIN_FOREACH{
#line 1016
foreach_visible (view)
 if (cfilter (point, d, cmin)) {
   coord n = facet_normal (point, d, fs);
   double alpha = plane_alpha (val(d,0,0,0), n);
   coord v[12];
   int m = facets (n, alpha, v, larger);
   if (m > 2) {
     vertex_buffer_glBegin (0x0009);
     for (int i = 0; i < m; i++) {
       if (linear) {
  if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, interp (point, v[i], col), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v[i], col); if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
       }
       else {
  if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
       }
       glnormal3d (view, n.x, n.y, n.z);
       glvertex3d (view,
     x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
     }
     vertex_buffer_glEnd ();
     view->ni++;
   }
 }end_foreach_visible();}END_FOREACH 
    }end_colorized();}
#line 1050 "/home/esteban/ProgramFile/basilisk/src/draw.h"
  if (expr) delete(((scalar[]){col,{-1}}));
  {end_tracing("draw_vof","/home/esteban/ProgramFile/basilisk/src/draw.h",1051);return true;}
end_tracing("draw_vof","/home/esteban/ProgramFile/basilisk/src/draw.h",1052);}
#line 1063 "/home/esteban/ProgramFile/basilisk/src/draw.h"
     
bool isoline (char * phi,
       double val,
       int n,
       bool edges,
       double larger, int filled,
       char * color,
       double min, double max, double spread,
       bool linear,
       Colormap map,
       float fc[3], float lc[3], float lw,
       bool expr,
       bool cbar, float size, float pos[2], char * label, double lscale, bool horizontal, bool border, bool mid, float fsize, char * format, int levels)
{tracing("isoline","/home/esteban/ProgramFile/basilisk/src/draw.h",1064);
#line 1104 "/home/esteban/ProgramFile/basilisk/src/draw.h"
  if (!(false)) qassert ("/home/esteban/ProgramFile/basilisk/src/draw.h", 1104, "false");

  {end_tracing("isoline","/home/esteban/ProgramFile/basilisk/src/draw.h",1106);return true;}
end_tracing("isoline","/home/esteban/ProgramFile/basilisk/src/draw.h",1107);}
#line 1120 "/home/esteban/ProgramFile/basilisk/src/draw.h"
     
bool cells (coord n, double alpha,
     float lc[3], float lw)
{tracing("cells","/home/esteban/ProgramFile/basilisk/src/draw.h",1121);
  bview * view = draw();
  {begin_draw_lines (view, lc, lw); {
#line 1137 "/home/esteban/ProgramFile/basilisk/src/draw.h"
     BEGIN_FOREACH{foreach_visible_plane (view, n, alpha) {
      coord v[12];
      int m = facets (n, alpha, v, 1.);
      if (m > 2) {
 vertex_buffer_glBegin (0x0002);
 for (int i = 0; i < m; i++)
   glvertex3d (view, x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
 vertex_buffer_glEnd ();
 view->ni++;
      }
    }end_foreach_visible_plane();}END_FOREACH 

  }end_draw_lines();}
  {end_tracing("cells","/home/esteban/ProgramFile/basilisk/src/draw.h",1150);return true;}
end_tracing("cells","/home/esteban/ProgramFile/basilisk/src/draw.h",1151);}






     
bool vectors (char * u, double scale, float lc[3], float lw)
{tracing("vectors","/home/esteban/ProgramFile/basilisk/src/draw.h",1159);
#line 1197 "/home/esteban/ProgramFile/basilisk/src/draw.h"
  fprintf (ferr, "vectors() is not implemented in 3D yet\n");

  {end_tracing("vectors","/home/esteban/ProgramFile/basilisk/src/draw.h",1199);return true;}
end_tracing("vectors","/home/esteban/ProgramFile/basilisk/src/draw.h",1200);}
#line 1220 "/home/esteban/ProgramFile/basilisk/src/draw.h"
     
bool squares (char * color,
       char * z,
       double min, double max, double spread,
       bool linear,
       Colormap map,
       float fc[3], float lc[3],
       bool expr,

       coord n,
       double alpha,
       float lw,
       bool cbar, float size, float pos[2], char * label, double lscale, bool horizontal, bool border, bool mid, float fsize, char * format, int levels)
{tracing("squares","/home/esteban/ProgramFile/basilisk/src/draw.h",1221);
#line 1248 "/home/esteban/ProgramFile/basilisk/src/draw.h"
  scalar col = {-1}; if (color && strcmp (color, "level")) { col = compile_expression (color, &expr); if (col.i < 0) {end_tracing("squares","/home/esteban/ProgramFile/basilisk/src/draw.h",1248);return false;} } double cmap[127][3]; if (color) { if (min == 0 && max == 0) { if (col.i < 0) min = 0, max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (spread < 0.) min = s.min, max = s.max; else { if (!spread) spread = 5.; min = avg - spread*s.stddev; max = avg + spread*s.stddev; } } } if (!map) map = jet; (* map) (cmap); } if ((3 > 2 || linear) && !fc[0] && !fc[1] && !fc[2]) fc[0] = fc[1] = fc[2] = 1.; if (cbar) colorbar (map, size, pos, label, lscale, min, max, horizontal, border, mid, lc, lw, fsize, format, levels);;
  scalar f = col;

  bview * view = draw();
  glShadeModel (0x1D01);
  if (linear) {
    {begin_colorized (fc, !VertexBuffer.color || !color, cmap, !VertexBuffer.color && color && linear && col.i >= 0); {
#line 1306 "/home/esteban/ProgramFile/basilisk/src/draw.h"
       BEGIN_FOREACH{foreach_visible_plane (view, n, alpha)
 if (val(f,0,0,0) != 1e30) {
   coord v[12];
   int m = facets (n, alpha, v, 1.);
   if (m > 2) {
     coord c = {0,0,0};
     for (int i = 0; i < m; i++)
       {
  c.x += v[i].x/m;
  
#line 1314
c.y += v[i].y/m;
  
#line 1314
c.z += v[i].z/m;}
     vertex_buffer_glBegin (0x0006);
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, interp (point, c, f), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, c, f); if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex3d (view, x + c.x*Delta, y + c.y*Delta, z + c.z*Delta);
     for (int i = 0; i < m; i++) {
       if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, interp (point, v[i], f), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v[i], f); if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
       glvertex3d (view,
     x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
     }
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, interp (point, v[0], f), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v[0], f); if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
     glvertex3d (view,
   x + v[0].x*Delta, y + v[0].y*Delta, z + v[0].z*Delta);
     vertex_buffer_glEnd ();
     view->ni++;
   }
 }end_foreach_visible_plane();}END_FOREACH 

    }end_colorized();}
  }
  else {
#line 1351 "/home/esteban/ProgramFile/basilisk/src/draw.h"
     BEGIN_FOREACH{foreach_visible_plane (view, n, alpha)
      if (val(f,0,0,0) != 1e30) {
 coord v[12];
 int m = facets (n, alpha, v, 1.);
 if (m > 2) {
   vertex_buffer_glBegin (0x0009);
   for (int i = 0; i < m; i++) {
     if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
     glvertex3d (view,
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
   }
   vertex_buffer_glEnd ();
   view->ni++;
 }
      }end_foreach_visible_plane();}END_FOREACH 

  }
  if (expr) delete (((scalar[]){col,{-1}}));




  {end_tracing("squares","/home/esteban/ProgramFile/basilisk/src/draw.h",1373);return true;}
end_tracing("squares","/home/esteban/ProgramFile/basilisk/src/draw.h",1374);}
#line 1385 "/home/esteban/ProgramFile/basilisk/src/draw.h"
     
bool box (bool notics, float lc[3], float lw)
{tracing("box","/home/esteban/ProgramFile/basilisk/src/draw.h",1386);
  bview * view = draw();
  {begin_draw_lines (view, lc, lw); {

    float height = 0.5*gl_StrokeHeight();
    float width = gl_StrokeWidth ('1'), scale = L0/(60.*width), length;
    float Z1 = 3 == 2 ? 0. : Z0;
    char label[80];

    glMatrixMode (0x1700);

    if (!notics) {
      int nt = 8;
      for (int i = 0; i <= nt; i++) {
 glPushMatrix();
 glTranslatef (X0 + i*L0/nt - height/2.*scale, Y0 - width/3.*scale, Z1);
 glRotatef (-90, 0, 0, 1);
 glScalef (scale, scale, 1.);
 sprintf (label, "%g", X0 + i*L0/nt);
 gl_StrokeString (label);
 glPopMatrix();

 glPushMatrix();
 sprintf (label, "%g", Y0 + i*L0/nt);
 length = gl_StrokeLength (label);
 glTranslatef (X0 - (length + width/3.)*scale,
        Y0 + i*L0/nt - height/2.*scale, Z1);
 glScalef (scale, scale, 1.);
 gl_StrokeString (label);
 glPopMatrix();


 glPushMatrix();
 sprintf (label, "%g", Z0 + i*L0/nt);
 length = gl_StrokeLength (label);
 glTranslatef (X0 - (length + width/3.)*scale,
        Y0, Z0 + i*L0/nt + height/2.*scale);
 glRotatef (-90, 1, 0, 0);
 glScalef (scale, scale, 1.);
 gl_StrokeString (label);
 glPopMatrix();

      }

      glPushMatrix();
      sprintf (label, "%g", X0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 + L0/2 - height*scale, Y0 - (length + 4.*width)*scale, Z1);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("X");
      glPopMatrix();


      glPushMatrix();
      sprintf (label, "%g", Y0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 - (length + 4.*width)*scale,
      Y0 + L0/2. - height*scale, Z1);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("Y");
      glPopMatrix();


      glPushMatrix();
      sprintf (label, "%g", Z0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 - (length + 4.*width)*scale,
      Y0, Z0 + L0/2. + height*scale);
      glRotatef (-90, 1, 0, 0);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("Z");
      glPopMatrix();

    }
#line 1473 "/home/esteban/ProgramFile/basilisk/src/draw.h"
    
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
 BEGIN_FOREACH{
#line 1473
foreach_level (0) {
      for (int i = -1; i <= 1; i += 2) {
 vertex_buffer_glBegin (0x0002);
 glvertex3d (view, x - Delta_x/2., y - Delta_y/2., z + i*Delta/2.);
 glvertex3d (view, x + Delta_x/2., y - Delta_y/2., z + i*Delta/2.);
 glvertex3d (view, x + Delta_x/2., y + Delta_y/2., z + i*Delta/2.);
 glvertex3d (view, x - Delta_x/2., y + Delta_y/2., z + i*Delta/2.);
 vertex_buffer_glEnd ();
 view->ni++;
 vertex_buffer_glBegin (0x0001);
 for (int j = -1; j <= 1; j += 2) {
   glvertex3d (view, x + i*Delta/2., y + j*Delta/2., z - Delta/2.);
   glvertex3d (view, x + i*Delta/2., y + j*Delta/2., z + Delta/2.);
 }
 vertex_buffer_glEnd ();
 view->ni++;
      }
    }end_foreach_level();}END_FOREACH 
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif


  
#line 1492
}end_draw_lines();}
  {end_tracing("box","/home/esteban/ProgramFile/basilisk/src/draw.h",1493);return true;}
end_tracing("box","/home/esteban/ProgramFile/basilisk/src/draw.h",1494);}
#line 1507 "/home/esteban/ProgramFile/basilisk/src/draw.h"
     
bool isosurface (char * f,
   double v,

   char * color,
   double min, double max, double spread,
   bool linear,
   Colormap map,
   float fc[3], float lc[3], float lw,
   bool expr,
   bool cbar, float size, float pos[2], char * label, double lscale, bool horizontal, bool border, bool mid, float fsize, char * format, int levels)
{tracing("isosurface","/home/esteban/ProgramFile/basilisk/src/draw.h",1508);

  if (!f)
    {end_tracing("isosurface","/home/esteban/ProgramFile/basilisk/src/draw.h",1521);return false;}

  scalar ff = {-1};
  bool fexpr;
  if (strcmp (f, "level")) {
    ff = compile_expression (f, &fexpr);
    if (ff.i < 0)
      {end_tracing("isosurface","/home/esteban/ProgramFile/basilisk/src/draw.h",1528);return false;}
  }

  scalar col = {-1}; if (color && strcmp (color, "level")) { col = compile_expression (color, &expr); if (col.i < 0) {end_tracing("isosurface","/home/esteban/ProgramFile/basilisk/src/draw.h",1531);return false;} } double cmap[127][3]; if (color) { if (min == 0 && max == 0) { if (col.i < 0) min = 0, max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (spread < 0.) min = s.min, max = s.max; else { if (!spread) spread = 5.; min = avg - spread*s.stddev; max = avg + spread*s.stddev; } } } if (!map) map = jet; (* map) (cmap); } if ((3 > 2 || linear) && !fc[0] && !fc[1] && !fc[2]) fc[0] = fc[1] = fc[2] = 1.; if (cbar) colorbar (map, size, pos, label, lscale, min, max, horizontal, border, mid, lc, lw, fsize, format, levels);;

  scalar  fv=new_vertex_scalar("fv");
  foreach_vertex_stencil(1,{0},)
    {_stencil_val(ff,0,0,0); _stencil_val(ff,-1,0,0); _stencil_val(ff,0,-1,0); _stencil_val(ff,-1,-1,0);
     _stencil_val(ff,0,0,-1); _stencil_val(ff,-1,0,-1); _stencil_val(ff,0,-1,-1); _stencil_val(ff,-1,-1,-1);
#line 1535
_stencil_val_a(fv,0,0,0);         
}end_foreach_vertex_stencil()
   BEGIN_FOREACH{
#line 1534
foreach_vertex()
    val(fv,0,0,0) = (val(ff,0,0,0) + val(ff,-1,0,0) + val(ff,0,-1,0) + val(ff,-1,-1,0) +
     val(ff,0,0,-1) + val(ff,-1,0,-1) + val(ff,0,-1,-1) + val(ff,-1,-1,-1))/8.;end_foreach_vertex();}END_FOREACH 

  vector  n=new_vector("n");
  foreach_stencil(1,{0},)
    {
      {_stencil_val(ff,1,0,0); _stencil_val(ff,-1,0,0);_stencil_val_a(n.x,0,0,0);   }
      
#line 1541
{_stencil_val(ff,0,1,0); _stencil_val(ff,0,-1,0);_stencil_val_a(n.y,0,0,0);   }
      
#line 1541
{_stencil_val(ff,0,0,1); _stencil_val(ff,0,0,-1);_stencil_val_a(n.z,0,0,0);   }}end_foreach_stencil()
   BEGIN_FOREACH{
#line 1539
foreach()
    {
      val(n.x,0,0,0) = ((val(ff,1,0,0) - val(ff,-1,0,0))/(2.*Delta));
      
#line 1541
val(n.y,0,0,0) = ((val(ff,0,1,0) - val(ff,0,-1,0))/(2.*Delta));
      
#line 1541
val(n.z,0,0,0) = ((val(ff,0,0,1) - val(ff,0,0,-1))/(2.*Delta));}end_foreach();}END_FOREACH 

  bview * view = draw();
  glShadeModel (0x1D01);
  {begin_colorized (fc, !VertexBuffer.color || !color, cmap, !VertexBuffer.color && color && linear && col.i >= 0); {
    foreach_visible_stencil (1,DEPAREN({view}),) {  
       
       
_stencil_val(fv,0,1,1); _stencil_val(fv,1,1,1); _stencil_val(fv,1,1,0);
 
#line 1549
_stencil_val(fv,0,1,0); 
#line 1548
_stencil_val(fv,0,0,1); _stencil_val(fv,1,0,1); _stencil_val(fv,1,0,0);
 
#line 1548
_stencil_val(fv,0,0,0);       
      
         
            

{
 if (color && (!linear || col.i < 0)) {               _stencil_val(col,0,0,0);     } 
 
 for (int j = 0; j < 3; j++) {      
   
   
     { _stencil_interp (point,NULL , n.x); }
     
#line 1559
{ _stencil_interp (point,NULL , n.y); }
     
#line 1559
{ _stencil_interp (point,NULL , n.z); }    
   
   if (linear) {
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) {        _stencil_interp (point,NULL , col);     } else {    _stencil_interp (point,NULL , col);                } }
   }
   else {
     if (color && (!linear || col.i < 0)) {               _stencil_val(col,0,0,0);     }
   }          
   
 } 
 
 
      }
    }end_foreach_visible_stencil()
     BEGIN_FOREACH{
#line 1546
foreach_visible (view) {
      double val[8] = {
 val(fv,0,0,0), val(fv,1,0,0), val(fv,1,0,1), val(fv,0,0,1),
 val(fv,0,1,0), val(fv,1,1,0), val(fv,1,1,1), val(fv,0,1,1)
      };
      double t[5][3][3];
      int nt = polygonize (val, v, t);
      for (int i = 0; i < nt; i++) {
 if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 vertex_buffer_glBegin (0x0009);
 for (int j = 0; j < 3; j++) {
   coord v = {t[i][j][0], t[i][j][1], t[i][j][2]}, np;
   
     np.x = interp (point, v, n.x);
     
#line 1559
np.y = interp (point, v, n.y);
     
#line 1559
np.z = interp (point, v, n.z);
   glnormal3d (view, np.x, np.y, np.z);
   if (linear) {
     if (color && linear && col.i >= 0) { if (VertexBuffer.color) { Color b = colormap_color (cmap, interp (point, v, col), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v, col); if (max > min) glTexCoord1d (clamp(((_v) - min)/(max - min), 0., 1.)); else glTexCoord1d (0.); } };
   }
   else {
     if (color && (!linear || col.i < 0)) { Color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), min, max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
   }
   glvertex3d (view, x + v.x*Delta_x, y + v.y*Delta_y, z + v.z*Delta_z);
 }
 vertex_buffer_glEnd ();
 view->ni++;
      }
    }end_foreach_visible();}END_FOREACH 
  }end_colorized();}
  if (expr) delete (((scalar[]){col,{-1}}));
  if (fexpr) delete (((scalar[]){ff,{-1}}));

  {delete((scalar*)((scalar[]){n.x,n.y,n.z,fv,{-1}}));{end_tracing("isosurface","/home/esteban/ProgramFile/basilisk/src/draw.h",1577);return true;}}delete((scalar*)((scalar[]){n.x,n.y,n.z,fv,{-1}}));
end_tracing("isosurface","/home/esteban/ProgramFile/basilisk/src/draw.h",1578);}
#line 1591 "/home/esteban/ProgramFile/basilisk/src/draw.h"
void travelling (double start, double end,
   float tx, float ty, float quat[4], float fov)
{
  static float stx, sty, squat[4], sfov;
  static double told = -1.;
  if (told < start && t >= start) {
    bview * view = get_view();
    stx = view->tx, sty = view->ty, sfov = view->fov;
    for (int i = 0; i < 4; i++)
      squat[i] = view->quat[i];
  }
  if (t >= start && t <= end)
    view ( (!tx ? stx : ((t - start)*(tx) + (end - t)*(stx))/(end - start)), (!ty ? sty : ((t - start)*(ty) + (end - t)*(sty))/(end - start))
, (!fov ? sfov : ((t - start)*(fov) + (end - t)*(sfov))/(end - start))
,
#line 51
(
    
#line 51
float[4]) 
#line 1605
{(!quat[0] ? squat[0] : ((t - start)*(quat[0]) + (end - t)*(squat[0]))/(end - start)), (!quat[1] ? squat[1] : ((t - start)*(quat[1]) + (end - t)*(squat[1]))/(end - start)),
           (!quat[2] ? squat[2] : ((t - start)*(quat[2]) + (end - t)*(squat[2]))/(end - start)), (!quat[3] ? squat[3] : ((t - start)*(quat[3]) + (end - t)*(squat[3]))/(end - start))}
#line 51
, 
1., 1., 1., 
800, 800, 4,
(
    
#line 54
float[3]) {0}, 
0., 0., 0., 
false, 
0., 0., 0., 
0., 
NULL, 
NULL, 
0, 
0., 0., 0., 0., 
NULL
#line 1606
);
  if (told < end && t >= end) {
    bview * view = get_view();
    stx = view->tx, sty = view->ty, sfov = view->fov;
    for (int i = 0; i < 4; i++)
      squat[i] = view->quat[i];
  }
  told = t;
}
#line 1631 "/home/esteban/ProgramFile/basilisk/src/draw.h"
     
bool draw_string (char * str,
    int pos,
    float size,
    float lc[3], float lw)

{tracing("draw_string","/home/esteban/ProgramFile/basilisk/src/draw.h",1632);
  bview * view = draw();

  glMatrixMode (0x1701);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode (0x1700);
  glPushMatrix();
  glLoadIdentity();

  vertex_buffer_glColor3f (lc[0], lc[1], lc[2]);
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));

  float width = gl_StrokeWidth ('1'), height = gl_StrokeHeight();
  if (!size)
    size = 40;
  float hscale = 2./(size*width), vscale = hscale*view->width/view->height;
  float vmargin = width/2.*vscale;
  if (pos == 0)
    glTranslatef (-1., -1. + vmargin, 0.);
  else if (pos == 1)
    glTranslatef (-1., 1. - height*vscale, 0.);
  else if (pos == 2)
    glTranslatef (1. - strlen(str)*width*hscale, 1. - height*vscale, 0.);
  else
    glTranslatef (1. - strlen(str)*width*hscale, -1. + vmargin, 0.);
  glScalef (hscale, vscale, 1.);
  gl_StrokeString (str);

  glMatrixMode (0x1700);
  glPopMatrix();
  glMatrixMode (0x1701);
  glPopMatrix();

  {end_tracing("draw_string","/home/esteban/ProgramFile/basilisk/src/draw.h",1672);return true;}
end_tracing("draw_string","/home/esteban/ProgramFile/basilisk/src/draw.h",1673);}




     
bool labels (char * f,
      float lc[3], float lw)
{tracing("labels","/home/esteban/ProgramFile/basilisk/src/draw.h",1679);
#line 1710 "/home/esteban/ProgramFile/basilisk/src/draw.h"
  fprintf (ferr, "labels() is not implemented in 3D yet\n");
  {end_tracing("labels","/home/esteban/ProgramFile/basilisk/src/draw.h",1711);return false;}

end_tracing("labels","/home/esteban/ProgramFile/basilisk/src/draw.h",1713);}







#line 1 "draw_json.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/draw_json.h"

int _view_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"view\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"tx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"ty\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fov\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"quat\": { \"type\": \"pfloat\", \"cardinality\": 4, \"value\": [%f,%f,%f,%f] }", 0., 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sy\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sz\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"width\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"height\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"samples\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"bg\": { \"type\": \"pfloat\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"theta\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"phi\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"psi\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"relative\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"tz\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"near\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"far\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"res\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"camera\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"cache\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p1x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p1y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p2x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p2y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _begin_translate_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"begin_translate\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"z\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _begin_mirror_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"begin_mirror\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"n\": { \"type\": \"pdouble\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _draw_vof_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"draw_vof\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"c\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"s\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"edges\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"larger\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"filled\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _isoline_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"isoline\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"phi\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"val\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"n\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"edges\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"larger\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"filled\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _cells_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"cells\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"n\": { \"type\": \"pdouble\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _vectors_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"vectors\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"u\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"scale\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _squares_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"squares\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"z\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"n\": { \"type\": \"pdouble\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _box_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"box\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"notics\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _isosurface_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"isosurface\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"f\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"v\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _travelling_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"travelling\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"start\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"end\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"tx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"ty\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"quat\": { \"type\": \"pfloat\", \"cardinality\": 4, \"value\": [%f,%f,%f,%f] }", 0., 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fov\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _draw_string_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"draw_string\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"str\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"pos\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", 0);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"size\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
int _labels_json (char * s, int len) {
  int i, len1 = 0;
  i = snprintf (s, len, "  \"labels\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"f\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", "");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", 0., 0., 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", 0.);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
#line 1722 "/home/esteban/ProgramFile/basilisk/src/draw.h"

struct {
  int (* json) (char * s, int len);
} bview_interface[] = {
  { _draw_vof_json },
  { _squares_json },
  { _cells_json },
  { _box_json },





  { _isosurface_json },

  { NULL }
};
#line 422 "/home/esteban/ProgramFile/basilisk/src/view.h"
#line 437 "/home/esteban/ProgramFile/basilisk/src/view.h"
bool load (FILE * fp, char * file, Array * buf);

static void bview_draw (bview * view)
{
  if (!view->active)
    return;
  view->active = false;
  glFinish ();
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#line 494 "/home/esteban/ProgramFile/basilisk/src/view.h"
     
bool save (char * file, char * format, char * opt,
    FILE * fp,
    float lw,
    int sort, int options,

    bview * view)
{tracing("save","/home/esteban/ProgramFile/basilisk/src/view.h",495);
  if (file) {
    char * s = strchr (file, '.'), * dot = s;
    while (s) {
      dot = s;
      s = strchr (s + 1, '.');
    }
    if (dot)
      format = dot + 1;
  }

  if (!view)
    view = get_view();

  if ((!strcmp (format, "png") && which ("convert")) ||
      !strcmp (format, "jpg") ||
      (file && is_animation (file))) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0) {
      FILE * fp = open_image (file, opt);
      if (!fp) {
 perror (file);
 {end_tracing("save","/home/esteban/ProgramFile/basilisk/src/view.h",524);return false;}
      }
      gl_write_image (fp, image, view->width, view->height, view->samples);
      close_image (file, fp);
    }
    {end_tracing("save","/home/esteban/ProgramFile/basilisk/src/view.h",529);return true;}
  }

  if (file && (fp = fopen (file, "w")) == NULL) {
    perror (file);
    {end_tracing("save","/home/esteban/ProgramFile/basilisk/src/view.h",534);return false;}
  }
  if (!fp)
    fp = fout;

  if (!strcmp (format, "ppm")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image (fp, image, view->width, view->height, view->samples);
  }

  else if (!strcmp (format, "png")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image_png (fp, image, view->width, view->height, view->samples);
  }

  else if (!strcmp (format, "bv")) {

    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
      format);
#line 573 "/home/esteban/ProgramFile/basilisk/src/view.h"
  }

  else if (!strcmp (format, "gnu") ||
    !strcmp (format, "obj") ||
    !strcmp (format, "kml") ||
    !strcmp (format, "ps") ||
    !strcmp (format, "eps") ||
    !strcmp (format, "tex") ||
    !strcmp (format, "pdf") ||
    !strcmp (format, "svg") ||
    !strcmp (format, "pgf"))
    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
      format);

  else {
    fprintf (ferr, "save(): unknown format '%s'\n", format);
    if (file) {
      fclose (fp);
      remove (file);
    }
    {end_tracing("save","/home/esteban/ProgramFile/basilisk/src/view.h",593);return false;}
  }

  fflush (fp);
  if (file)
    fclose (fp);

  {end_tracing("save","/home/esteban/ProgramFile/basilisk/src/view.h",600);return true;}
end_tracing("save","/home/esteban/ProgramFile/basilisk/src/view.h",601);}







static char * remove_blanks (char * line)
{
  while (strchr (" \t", *line)) line++;
  char * s = line, * cur = line;
  bool instring = false;
  while (*s != '\0' && *s != '#') {
    if (*s == '"')
      instring = !instring;
    if (instring || !strchr (" \t", *s))
      *cur++ = *s;
    s++;
  }
  *cur = '\0';
  return line;
}

#line 1 "parse.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/parse.h"




enum ParamsType { pstring, pint, punsigned, pbool, pfloat, pdouble, pcolormap };

typedef struct {
  char * key;
  enum ParamsType type;
  void * val;
  int n;
} Params;

static bool atobool (char * s)
{
  if (!strcmp (s, "true"))
    return true;
  if (!strcmp (s, "false"))
    return false;
  return atoi (s) != 0;
}

static bool args (Params * p, char * val)
{
  static char * name[] = { "string", "int", "unsigned",
      "bool", "float", "double", "colormap" };
  switch (p->type) {

  case pstring:
    if (val[0] != '"') {
      fprintf (ferr, "expecting a string for '%s' got '%s'\n", p->key, val);
      return false;
    }
    if (val[strlen(val) - 1] != '"') {
      fprintf (ferr, "unterminated quoted string '%s'\n", val);
      return false;
    }
    val[strlen(val) - 1] = '\0';
    char * s = &val[1];
    int nc = 0;
    while (*s != '\0') {
      if (!strchr (" \t\n\r", *s))
 nc++;
      s++;
    }
    *((char **)p->val) = nc > 0 ? &val[1] : NULL;
    break;

  case pcolormap:
    if (!strcmp (val, "jet"))
      *((Colormap *)p->val) = jet;
    else if (!strcmp (val, "cool_warm"))
      *((Colormap *)p->val) = cool_warm;
    else if (!strcmp (val, "gray"))
      *((Colormap *)p->val) = gray;
    else if (!strcmp (val, "randomap"))
      *((Colormap *)p->val) = randomap;
    else {
      fprintf (ferr, "unknown colormap '%s'\n", val);
      return false;
    }
    break;

  case pint: case punsigned: case pbool: case pdouble: case pfloat:
    if (val[0] == '"') {
      fprintf (ferr, "expecting a %s for '%s' got %s\n",
        name[p->type], p->key, val);
      return false;
    }
    if (!p->n) {
      switch (p->type) {
      case pint: *((int *)p->val) = atoi(val); break;
      case punsigned: *((unsigned *)p->val) = atoi(val); break;
      case pbool: *((bool *)p->val) = atobool(val); break;
      case pfloat: *((float *)p->val) = atof(val); break;
      case pdouble: *((double *)p->val) = atof(val); break;
      default: if (!(false)) qassert ("/home/esteban/ProgramFile/basilisk/src/parse.h", 77, "false");
      }
    }
    else {
      if (val[0] != '{') {
 fprintf (ferr, "expecting an array for '%s' got %s\n", p->key, val);
 return false;
      }
      val++;
      int i = 0;
      char c = ',';
      while (i < p->n && c != '}') {
 char * s = strchr (val, ',');
 if (!s)
   s = strchr (val, '}');
 if (!s) {
   fprintf (ferr, "expecting an array for '%s' got %s\n", p->key, val);
   return false;
 }
 c = *s;
 *s++ = '\0';
 switch (p->type) {
 case pint: ((int *)p->val)[i++] = atoi (val); break;
 case punsigned: ((unsigned *)p->val)[i++] = atoi (val); break;
 case pbool: ((bool *)p->val)[i++] = atobool (val); break;
 case pfloat: ((float *)p->val)[i++] = atof (val); break;
 case pdouble: ((double *)p->val)[i++] = atof (val); break;
 default: if (!(false)) qassert ("/home/esteban/ProgramFile/basilisk/src/parse.h", 104, "false");
 }
 val = s;
      }
      if (c != '}') {
 fprintf (ferr, "expecting '}' for '%s' got %s\n", p->key, val);
 return false;
      }
    }
    break;

  default:
    if (!(false)) qassert ("/home/esteban/ProgramFile/basilisk/src/parse.h", 116, "false");
  }
  return true;
}

static char * find_comma (char * s)
{
  int par = 0;
  while (*s != '\0') {
    if (*s == ',' && par == 0) {
      *s = '\0';
      return s + 1;
    }
    if (*s == '{')
      par++;
    else if (*s == '}')
      par--;
    s++;
  }
  return NULL;
}

static char * mystrtok (char * str, const char * delim)
{
  static char * s = NULL;
  char * start = str ? str : s;
  bool string = false;
  s = start;
  while (*s != '\0') {
    if (*s == '"')
      string = !string;
    if (!string && strchr(delim, *s))
      break;
    s++;
  }
  if (*s != '\0')
    *s++ = '\0';
  return start;
}

int parse_params (Params * params)
{
  char * s;
  int i = 0, n = 0;
  Params * p = params;
  while (p->key) p++, n++;
  if (!(s = mystrtok (NULL, ");")) || s[0] == '\n')
    return false;
  while (s && *s != '\0') {
    char * next = find_comma (s), * key = s;
    if ((s = strchr (key, '='))) {
      s[0] = '\0', s++;
      i = -1;
      Params * p = params;
      while (p->key && strcmp(p->key, key)) p++;
      if (!p->key) {
 fprintf (ferr, "unknown key '%s'\n", key);
 return false;
      }
      if (!args (p, s))
 return false;
    }
    else {
      if (i < 0) {
 fprintf (ferr, "anonymous value '%s' after keys\n", key);
 return false;
      }
      if (i >= n) {
 fprintf (ferr, "too many parameters: '%s' %d %d\n", key, i, n);
 return false;
      }
      if (!args (&params[i], key))
 return false;
      i++;
    }
    s = next;
  }
  return true;
}
#line 626 "/home/esteban/ProgramFile/basilisk/src/view.h"

bool process_line (char * line)
{
  if (line[0] == '\0')
    return true;
  char * s = mystrtok (remove_blanks (line), "(");
  if (!s)
    return true;

  if (!strcmp (s, "restore")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      bview * view = get_view();
      if (view->cache) {
 free_cexpr (view->cache);
 view->cache = pcalloc (1, sizeof (cexpr),__func__,__FILE__,__LINE__);
      }
      if (!restore ( file, all
#line 1125 "/home/esteban/ProgramFile/basilisk/src/output.h"
, 
NULL
#line 644 "/home/esteban/ProgramFile/basilisk/src/view.h"
))
 fprintf (ferr, "could not restore from '%s'\n", file);
      else {
 restriction (all);
 fields_stats();
 clear();
      }
    }
  }

  else if (!strcmp (s, "dump")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    dump ( file
#line 1047 "/home/esteban/ProgramFile/basilisk/src/output.h"
, 
all, 
NULL, 
false
#line 657 "/home/esteban/ProgramFile/basilisk/src/view.h"
);
  }

  else if (!strcmp (s, "input_gfs")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      input_gfs ( 
#line 157 "/home/esteban/ProgramFile/basilisk/src/input.h"
stdin
#line 664 "/home/esteban/ProgramFile/basilisk/src/view.h"
, all, file);
      restriction (all);
      fields_stats();
      clear();
    }
  }

  else if (!strcmp (s, "save")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      save ( file
#line 495
, "ppm", NULL, 
NULL, 
0, 
0, 0, 

NULL
#line 675
);
  }

  else if (!strcmp (s, "load")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      load ( 
#line 437
NULL
#line 682
, file
#line 437
, NULL
#line 682
);
  }






#line 1 "draw_get.h"
#line 1 "/home/esteban/ProgramFile/basilisk/src/draw_get.h"

else if (!strcmp (s, "view")) {
  float tx = 0.;
  float ty = 0.;
  float fov = 0.;
  float quat[4] = {0};
  float sx = 1.;
  float sy = 1.;
  float sz = 1.;
  unsigned width = 800;
  unsigned height = 800;
  unsigned samples = 4;
  float bg[3] = {0};
  float theta = 0.;
  float phi = 0.;
  float psi = 0.;
  bool relative = false;
  float tz = 0.;
  float near = 0.;
  float far = 0.;
  float res = 0.;
  char * camera = NULL;
  MapFunc map = NULL;
  int cache = 0;
  float p1x = 0.;
  float p1y = 0.;
  float p2x = 0.;
  float p2y = 0.;
  bview * view1 = NULL;
  Params params[] = {
    {"tx", pfloat, &tx},
    {"ty", pfloat, &ty},
    {"fov", pfloat, &fov},
    {"quat", pfloat, quat, 4},
    {"sx", pfloat, &sx},
    {"sy", pfloat, &sy},
    {"sz", pfloat, &sz},
    {"width", punsigned, &width},
    {"height", punsigned, &height},
    {"samples", punsigned, &samples},
    {"bg", pfloat, bg, 3},
    {"theta", pfloat, &theta},
    {"phi", pfloat, &phi},
    {"psi", pfloat, &psi},
    {"relative", pbool, &relative},
    {"tz", pfloat, &tz},
    {"near", pfloat, &near},
    {"far", pfloat, &far},
    {"res", pfloat, &res},
    {"camera", pstring, &camera},
    {"cache", pint, &cache},
    {"p1x", pfloat, &p1x},
    {"p1y", pfloat, &p1y},
    {"p2x", pfloat, &p2x},
    {"p2y", pfloat, &p2y},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  view (tx,ty,fov,quat,sx,sy,sz,width,height,samples,bg,theta,phi,psi,relative,tz,near,far,res,camera,map,cache,p1x,p1y,p2x,p2y,view1);
}
else if (!strcmp (s, "begin_translate")) {
  float x = 0;
  float y = 0.;
  float z = 0.;
  Params params[] = {
    {"x", pfloat, &x},
    {"y", pfloat, &y},
    {"z", pfloat, &z},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  begin_translate (x,y,z);
}
else if (!strcmp (s, "begin_mirror")) {
  coord n = {0};
  double alpha = 0.;
  Params params[] = {
    {"n", pdouble, &n, 3},
    {"alpha", pdouble, &alpha},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  begin_mirror (n,alpha);
}
else if (!strcmp (s, "draw_vof")) {
  char * c = NULL;
  char * s = NULL;
  bool edges = false;
  double larger = 0.;
  int filled = 0;
  char * color = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  float lw = 1.;
  bool expr = false;
  Params params[] = {
    {"c", pstring, &c},
    {"s", pstring, &s},
    {"edges", pbool, &edges},
    {"larger", pdouble, &larger},
    {"filled", pint, &filled},
    {"color", pstring, &color},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!draw_vof (c,s,edges,larger,filled,color,min,max,spread,linear,map,fc,lc,lw,expr
#line 893 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 122 "/home/esteban/ProgramFile/basilisk/src/draw_get.h"
))
    return false;
}
else if (!strcmp (s, "isoline")) {
  char * phi = NULL;
  double val = 0.;
  int n = 1;
  bool edges = false;
  double larger = 0.;
  int filled = 0;
  char * color = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  float lw = 1.;
  bool expr = false;
  Params params[] = {
    {"phi", pstring, &phi},
    {"val", pdouble, &val},
    {"n", pint, &n},
    {"edges", pbool, &edges},
    {"larger", pdouble, &larger},
    {"filled", pint, &filled},
    {"color", pstring, &color},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!isoline (phi,val,n,edges,larger,filled,color,min,max,spread,linear,map,fc,lc,lw,expr
#line 1074 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 161 "/home/esteban/ProgramFile/basilisk/src/draw_get.h"
))
    return false;
}
else if (!strcmp (s, "cells")) {
  coord n = {0,0,1};
  double alpha = 0.;
  float lc[3] = {0};
  float lw = 1.;
  Params params[] = {
    {"n", pdouble, &n, 3},
    {"alpha", pdouble, &alpha},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!cells (n,alpha,lc,lw))
    return false;
}
else if (!strcmp (s, "vectors")) {
  char * u = NULL;
  double scale = 1;
  float lc[3] = {0};
  float lw = 1.;
  Params params[] = {
    {"u", pstring, &u},
    {"scale", pdouble, &scale},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!vectors (u,scale,lc,lw))
    return false;
}
else if (!strcmp (s, "squares")) {
  char * color = NULL;
  char * z = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  bool expr = false;
  coord n = {0,0,1};
  double alpha = 0;
  Params params[] = {
    {"color", pstring, &color},
    {"z", pstring, &z},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"n", pdouble, &n, 3},
    {"alpha", pdouble, &alpha},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!squares (color,z,min,max,spread,linear,map,fc,lc,expr,n,alpha
#line 1230 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
1, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 226 "/home/esteban/ProgramFile/basilisk/src/draw_get.h"
))
    return false;
}
else if (!strcmp (s, "box")) {
  bool notics = false;
  float lc[3] = {0};
  float lw = 1.;
  Params params[] = {
    {"notics", pbool, &notics},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!box (notics,lc,lw))
    return false;
}
else if (!strcmp (s, "isosurface")) {
  char * f = NULL;
  double v = 0.;
  char * color = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  float lw = 1;
  bool expr = false;
  Params params[] = {
    {"f", pstring, &f},
    {"v", pdouble, &v},
    {"color", pstring, &color},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!isosurface (f,v,color,min,max,spread,linear,map,fc,lc,lw,expr
#line 1516 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 272 "/home/esteban/ProgramFile/basilisk/src/draw_get.h"
))
    return false;
}
else if (!strcmp (s, "travelling")) {
  double start = 0;
  double end = 0;
  float tx = 0;
  float ty = 0;
  float quat[4] = {0};
  float fov = 0;
  Params params[] = {
    {"start", pdouble, &start},
    {"end", pdouble, &end},
    {"tx", pfloat, &tx},
    {"ty", pfloat, &ty},
    {"quat", pfloat, quat, 4},
    {"fov", pfloat, &fov},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  travelling (start,end,tx,ty,quat,fov);
}
else if (!strcmp (s, "draw_string")) {
  char * str = NULL;
  int pos = 0;
  float size = 40;
  float lc[3] = {0};
  float lw = 1;
  Params params[] = {
    {"str", pstring, &str},
    {"pos", pint, &pos},
    {"size", pfloat, &size},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!draw_string (str,pos,size,lc,lw))
    return false;
}
else if (!strcmp (s, "labels")) {
  char * f = NULL;
  float lc[3] = {0};
  float lw = 1;
  Params params[] = {
    {"f", pstring, &f},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!labels (f,lc,lw))
    return false;
}
#line 691 "/home/esteban/ProgramFile/basilisk/src/view.h"

  else if (!strcmp (s, "end_mirror"))
    end_mirror();

  else if (!strcmp (s, "end_translate"))
    end_translate();

  else if (!strcmp (s, "clear"))
    clear();

  else if (s[0] != '\n' && s[0] != '\0')
    fprintf (ferr, "load(): syntax error: '%s'\n", s);

  return true;
}

bool load (FILE * fp, char * file, Array * buf)
{
  if (file) {
    fp = fopen (file, "r");
    if (!fp) {
      perror (file);
      return false;
    }
  }

  if (fp) {
    char line[256];
    while (fgets (line, 256, fp) && process_line (line));
  }
  else if (buf) {
    int i = 0;
    char * s = (char *) buf->p;
    while (i < buf->len) {
      char * start = s;
      while (i < buf->len && *s != '\n')
 s++, i++;
      if (*s == '\n' && ++s > start) {
 char line[s - start + 1];
 strncpy (line, start, s - start);
 line[s - start] = '\0';
 process_line (line);
      }
    }
  }
  return true;
}
#line 20 "main3D.c"
#line 28 "main3D.c"
const double tEnd = 0.1e0;
const double tStep = 5e-3;




const double jetThickness = 3e-3;
const double domainLength = 50. * jetThickness;





scalar  f0={16},  U0={17};
char name[80];


const double rhoL = 1000., rhoG = 1.225;

const double muL= 1e-3, muG = 1.8e-5;



const double gravity = 9.81;

const double sigma = 72e-3;


double Re;
const double We = 10;

double u0;
double lambda;




const int maxlevel = 6;
const double uemax = 1e-3;







static double _boundary6(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet(val(f0,0,0,0) * u0 * (1 - sq(y/jetThickness)), point, neighbor, _s, data));}}}static double _boundary6_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_homogeneous(val(f0,0,0,0) * u0 * (1 - sq(y/jetThickness)), point, neighbor, _s, data));}}}
static double _boundary7(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet(0., point, neighbor, _s, data));}}}static double _boundary7_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0., point, neighbor, _s, data));}}}
static double _boundary8(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet(0, point, neighbor, _s, data));}}}static double _boundary8_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0, point, neighbor, _s, data));}}}

static double _boundary9(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_face(val(f0,0,0,0) * u0 * (1 - sq(y/jetThickness)), point, neighbor, _s, data));}}}static double _boundary9_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_face_homogeneous(val(f0,0,0,0) * u0 * (1 - sq(y/jetThickness)), point, neighbor, _s, data));}}}
static double _boundary10(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet(0., point, neighbor, _s, data));}}}static double _boundary10_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0., point, neighbor, _s, data));}}}
static double _boundary11(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet(0., point, neighbor, _s, data));}}}static double _boundary11_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0., point, neighbor, _s, data));}}}

static double _boundary12(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( val(f0,0,0,0));}}}




static double _boundary13(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(0., point, neighbor, _s, data));}}}static double _boundary13_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(0., point, neighbor, _s, data));}}}
static double _boundary14(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(0., point, neighbor, _s, data));}}}static double _boundary14_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(0., point, neighbor, _s, data));}}}
static double _boundary15(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann(0., point, neighbor, _s, data));}}}static double _boundary15_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _neumann_homogeneous(0., point, neighbor, _s, data));}}}

static double _boundary16(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet(0., point, neighbor, _s, data));}}}static double _boundary16_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0., point, neighbor, _s, data));}}}
static double _boundary17(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet(0., point, neighbor, _s, data));}}}static double _boundary17_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0., point, neighbor, _s, data));}}}

static double _boundary18(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet(0., point, neighbor, _s, data));}}}static double _boundary18_homogeneous(Point point,Point neighbor,scalar _s,bool *data){int ig=0;NOT_UNUSED(ig);int jg=0;NOT_UNUSED(jg);int kg=0;NOT_UNUSED(kg);POINT_VARIABLES;{int ig=neighbor.i-point.i;if(ig==0)ig=_attribute[_s.i].d.x;NOT_UNUSED(ig);int jg=neighbor.j-point.j;if(jg==0)jg=_attribute[_s.i].d.y;NOT_UNUSED(jg);int kg=neighbor.k-point.k;if(kg==0)kg=_attribute[_s.i].d.z;NOT_UNUSED(kg);POINT_VARIABLES;{return( _dirichlet_homogeneous(0., point, neighbor, _s, data));}}}




int main()
{
#line 123
_init_solver();
    
#line 101
periodic(front);

    init_grid(1 << maxlevel);

    origin(0., -domainLength/2.
#line 105 "/home/esteban/ProgramFile/basilisk/src/common.h"
, 0.
#line 105 "main3D.c"
);
    size(domainLength);

    rho1 = rhoL;
    rho2 = rhoG;

    mu1 = muL;
    mu2 = muG;

    _attribute[f.i].sigma = sigma;

    u0 = sqrt((We * sigma) / (rhoL * jetThickness));


    Re = rhoL * u0 * jetThickness / muL;
    lambda = sqrt(sigma / ((rhoL-rhoG)*gravity)) * 1e3;

    run();
free_solver();

#line 123
}




static int init_0_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t = 0)!=0;*ip=i;*tp=t;return ret;}





#line 128
      static int init_0(const int i,const double t,Event *_ev){tracing("init_0","main3D.c",128);
{
    TOLERANCE = 1e-5;

    G.x = gravity;

    do { scalar  phi=new_vertex_scalar("phi"); foreach_vertex_stencil(1,{0},) {_stencil_val_a(phi,0,0,0);    }end_foreach_vertex_stencil()  BEGIN_FOREACH{foreach_vertex() val(phi,0,0,0) = sq(jetThickness) - sq(y);end_foreach_vertex();}END_FOREACH  fractions (phi, f0
#line 121 "/home/esteban/ProgramFile/basilisk/src/fractions.h"
,
(
  
#line 122
vector) {0}, 0.
#line 134 "main3D.c"
);delete((scalar*)((scalar[]){phi,{-1}})); } while(0);

    restriction(((scalar[]){f0,{-1}}));

    foreach_stencil (1,{0},)
    { 
_stencil_val(f0,0,0,0);
        
#line 140
_stencil_val_a(f,0,0,0);     
        _stencil_val_a(U0,0,0,0); 
_stencil_val(U0,0,0,0); _stencil_val(f,0,0,0);      
        
#line 142
_stencil_val_a(u.x,0,0,0);  
    }end_foreach_stencil()

     BEGIN_FOREACH{
#line 138
foreach ()
    {
        val(f,0,0,0) = val(f0,0,0,0) * (x < 0.015);
        val(U0,0,0,0) = u0 * (1 - sq(y/jetThickness));
        val(u.x,0,0,0) = val(U0,0,0,0) * val(f,0,0,0);
    }end_foreach();}END_FOREACH 
}{end_tracing("init_0","main3D.c",144);return 0;}end_tracing("init_0","main3D.c",144);}




static int logfile_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(i++)!=0;*ip=i;*tp=t;return ret;}





#line 149
      static int logfile(const int i,const double t,Event *_ev){tracing("logfile","main3D.c",149);
{
    if (i == 0)
    {
        fprintf(ferr, "-----------------------------------------------------\n");
        fprintf(ferr, "Re %.1f | We %.2f | u0 %.3f | g %.2f | l %.2f | D %.3f | s %1.1e \n", Re, We, u0, gravity, lambda, jetThickness*1000, sigma);
        fprintf(ferr, "-----------------------------------------------------\n");
        fprintf(ferr, " t / t_max | dt | p.i | f.i | u.i | ncells | perf \n");
        fprintf(ferr, "-----------------------------------------------------\n");

    }

    fprintf(ferr, "%.5f / %.5f | %.2e |  %2d |  %2d |  %2d | %6ld | %.2f\n", t, tEnd, dt, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t);
}{end_tracing("logfile","main3D.c",162);return 0;}end_tracing("logfile","main3D.c",162);}

static int logEnd_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t=tEnd)!=0;*ip=i;*tp=t;return ret;}


#line 164
      static int logEnd(const int i,const double t,Event *_ev){tracing("logEnd","main3D.c",164);
{
    fprintf(ferr, "-----------------------------------------------------\n");
    fprintf(ferr, "Re %.1f | We %.2f | u0 %.3f | g %.2f | l %.2f | D %.3f | s %1.1e \n", Re, We, u0, gravity, lambda, jetThickness*1000, sigma);
    fprintf(ferr, "-----------------------------------------------------\n");
}{end_tracing("logEnd","main3D.c",169);return 0;}end_tracing("logEnd","main3D.c",169);}
#line 182
static int movie_expr0(int *ip,double *tp,Event *_ev){int i=*ip;double t=*tp;int ret=(t += tStep)!=0;*ip=i;*tp=t;return ret;}
#line 182 "main3D.c"
      static int movie(const int i,const double t,Event *_ev){tracing("movie","main3D.c",182);
{
    view( 
#line 49 "/home/esteban/ProgramFile/basilisk/src/draw.h"
0.
#line 184 "main3D.c"
, .4, 38
#line 50 "/home/esteban/ProgramFile/basilisk/src/draw.h"
,
(
    
#line 51
float[4]) {0}, 
1., 1., 1.
#line 184 "main3D.c"
, 750, 750, 4,
#line 54 "/home/esteban/ProgramFile/basilisk/src/draw.h"
(
    
#line 54
float[3]) 
#line 184 "main3D.c"
{1,1,1}
#line 54 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
0., 0., 0., 
false, 
0., 0., 0., 
0.
#line 184 "main3D.c"
, "iso"
#line 59 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
NULL, 
0, 
0., 0., 0., 0., 
NULL
#line 184 "main3D.c"
);

    clear();


    squares ("u.x"
#line 1221 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
NULL, 
0, 0, 0
#line 189 "main3D.c"
, true
#line 1224 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
jet,
(
       
#line 1226
float[3]) {0},( float[3]) {0}, 
false
#line 189 "main3D.c"
,
#line 1229 "/home/esteban/ProgramFile/basilisk/src/draw.h"
(

       
#line 1229
coord) 
#line 189 "main3D.c"
{1,0,0}
#line 1229 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
0, 
1, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 189 "main3D.c"
);


    scalar  omega=new_scalar("omega");
    vorticity (u, omega);
    squares("omega"
#line 1221 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
NULL, 
0, 0, 0
#line 194 "main3D.c"
, true
#line 1224 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
jet,
(
       
#line 1226
float[3]) {0},( float[3]) {0}, 
false
#line 194 "main3D.c"
,
#line 1229 "/home/esteban/ProgramFile/basilisk/src/draw.h"
(

       
#line 1229
coord) 
#line 194 "main3D.c"
{0,1,0}
#line 1229 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, 
0, 
1, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 194 "main3D.c"
);


    draw_vof("f"
#line 886 "/home/esteban/ProgramFile/basilisk/src/draw.h"
, NULL, false, 
0., 0, 
NULL, 
0, 0, 0, 
false, 
jet,
(
        
#line 892
float[3]) {0},( float[3]) {0}, 1., 
false, 
false, 15,( float[2]) {-.95, -.95}, "", 1, false, false, false, 50, "%g", 50
#line 197 "main3D.c"
);

    save("movie.mp4"
#line 495 "/home/esteban/ProgramFile/basilisk/src/view.h"
, "ppm", NULL, 
NULL, 
0, 
0, 0, 

NULL
#line 199 "main3D.c"
);delete((scalar*)((scalar[]){omega,{-1}}));
}{end_tracing("movie","main3D.c",200);return 0;}end_tracing("movie","main3D.c",200);}
#line 2 "ast/init_solver.h"

static void _init_solver (void)
{
  void init_solver();
datasize=18*sizeof(real);
  
#line 6
init_solver();
  {
#line 24
multigrid_methods();

    

    
#line 12
{
      
      
    
      {  
#line 42 "/home/esteban/ProgramFile/basilisk/src/run.h"
event_register((Event){0,1,defaults,{defaults_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/run.h",42,"defaults"});  
#line 126 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
event_register((Event){0,1,defaults_0,{defaults_0_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",126,"defaults"});  
#line 194
event_register((Event){0,1,default_display,{default_display_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",194,"default_display"});  








event_register((Event){0,1,init,{init_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",203,"init"});  
#line 127 "/home/esteban/ProgramFile/basilisk/src/vof.h"
event_register((Event){0,1,defaults_1,{defaults_1_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/vof.h",127,"defaults"});  
#line 10 "/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h"
event_register((Event){0,1,defaults_2,{defaults_2_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h",10,"defaults"});  
#line 30 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
event_register((Event){0,1,defaults_3,{defaults_3_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/iforce.h",30,"defaults"});  
#line 37 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h"
event_register((Event){0,1,defaults_4,{defaults_4_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",37,"defaults"});  
#line 128 "main3D.c"
event_register((Event){0,1,init_0,{init_0_expr0},((int *)0),((double *)0),"main3D.c",128,"init"});  
#line 149
event_register((Event){0,1,logfile,{logfile_expr0},((int *)0),((double *)0),"main3D.c",149,"logfile"});  
#line 164
event_register((Event){0,1,logEnd,{logEnd_expr0},((int *)0),((double *)0),"main3D.c",164,"logEnd"});  
#line 182
event_register((Event){0,1,movie,{movie_expr0},((int *)0),((double *)0),"main3D.c",182,"movie"});
	
	
	
      
#line 22 "ast/init_solver.h"
}
#line 383 "/home/esteban/ProgramFile/basilisk/src/common.h"
init_const_vector((vector){{_NVARMAX+4},{_NVARMAX+5},{_NVARMAX+6}},"zerof",(double[]){0.,0.,0.});
init_const_vector((vector){{_NVARMAX+7},{_NVARMAX+8},{_NVARMAX+9}},"unityf",(double[]){1.,1.,1.});
init_const_scalar((scalar){_NVARMAX+10},"unity", 1.);
init_const_scalar((scalar){_NVARMAX+11},"zeroc", 0.);



init_const_vector((vector){{_NVARMAX+12},{_NVARMAX+13},{_NVARMAX+14}},"fm",(double[]){1.,1.,1.});
init_const_scalar((scalar){_NVARMAX+15},"cm", 1.);  init_scalar((scalar){0},"p");  init_vector((vector){{1},{2},{3}},"u");  init_vector((vector){{4},{5},{6}},"g");  init_scalar((scalar){7},"pf");  init_face_vector((vector){{8},{9},{10}},"uf");  init_scalar((scalar){11},"f");  init_face_vector((vector){{12},{13},{14}},"alphav");  init_scalar((scalar){15},"rhov");  init_scalar((scalar){16},"f0");  init_scalar((scalar){17},"U0");
    
#line 23 "ast/init_solver.h"
}_attribute[p.i].dirty=1,_attribute[p.i].boundary[right]=_boundary0,_attribute[p.i].boundary_homogeneous[right]=_boundary0_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[left]=_boundary1,_attribute[p.i].boundary_homogeneous[left]=_boundary1_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[top]=_boundary2,_attribute[p.i].boundary_homogeneous[top]=_boundary2_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[bottom]=_boundary3,_attribute[p.i].boundary_homogeneous[bottom]=_boundary3_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[front]=_boundary4,_attribute[p.i].boundary_homogeneous[front]=_boundary4_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[back]=_boundary5,_attribute[p.i].boundary_homogeneous[back]=_boundary5_homogeneous;_attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[left]=_boundary6,_attribute[u.x.i].boundary_homogeneous[left]=_boundary6_homogeneous;_attribute[u.y.i].dirty=1,_attribute[u.y.i].boundary[left]=_boundary7,_attribute[u.y.i].boundary_homogeneous[left]=_boundary7_homogeneous;_attribute[u.z.i].dirty=1,_attribute[u.z.i].boundary[left]=_boundary8,_attribute[u.z.i].boundary_homogeneous[left]=_boundary8_homogeneous;_attribute[uf.x.i].dirty=1,_attribute[uf.x.i].boundary[left]=_boundary9,_attribute[uf.x.i].boundary_homogeneous[left]=_boundary9_homogeneous;_attribute[uf.y.i].dirty=1,_attribute[uf.y.i].boundary[left]=_boundary10,_attribute[uf.y.i].boundary_homogeneous[left]=_boundary10_homogeneous;_attribute[uf.z.i].dirty=1,_attribute[uf.z.i].boundary[left]=_boundary11,_attribute[uf.z.i].boundary_homogeneous[left]=_boundary11_homogeneous;_attribute[f.i].dirty=1,_attribute[f.i].boundary[left]=_boundary12,_attribute[f.i].boundary_homogeneous[left]=_boundary12;_attribute[u.x.i].dirty=1,_attribute[u.x.i].boundary[right]=_boundary13,_attribute[u.x.i].boundary_homogeneous[right]=_boundary13_homogeneous;_attribute[u.y.i].dirty=1,_attribute[u.y.i].boundary[right]=_boundary14,_attribute[u.y.i].boundary_homogeneous[right]=_boundary14_homogeneous;_attribute[u.z.i].dirty=1,_attribute[u.z.i].boundary[right]=_boundary15,_attribute[u.z.i].boundary_homogeneous[right]=_boundary15_homogeneous;_attribute[p.i].dirty=1,_attribute[p.i].boundary[right]=_boundary16,_attribute[p.i].boundary_homogeneous[right]=_boundary16_homogeneous;_attribute[pf.i].dirty=1,_attribute[pf.i].boundary[right]=_boundary17,_attribute[pf.i].boundary_homogeneous[right]=_boundary17_homogeneous;_attribute[f.i].dirty=1,_attribute[f.i].boundary[right]=_boundary18,_attribute[f.i].boundary_homogeneous[right]=_boundary18_homogeneous;  
#line 50 "/home/esteban/ProgramFile/basilisk/src/run.h"
event_register((Event){0,1,cleanup,{cleanup_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/run.h",50,"cleanup"});  
#line 229 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h"
event_register((Event){0,1,set_dtmax,{set_dtmax_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",229,"set_dtmax"});  

event_register((Event){0,1,stability,{stability_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",231,"stability"});  









event_register((Event){0,1,vof,{vof_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",241,"vof"});  
event_register((Event){0,1,tracer_advection,{tracer_advection_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",242,"tracer_advection"});  
event_register((Event){0,1,tracer_diffusion,{tracer_diffusion_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",243,"tracer_diffusion"});  






event_register((Event){0,1,properties,{properties_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",250,"properties"});  
#line 323
event_register((Event){0,1,advection_term,{advection_term_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",323,"advection_term"});  
#line 352
event_register((Event){0,1,viscous_term,{viscous_term_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",352,"viscous_term"});  
#line 388
event_register((Event){0,1,acceleration,{acceleration_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",388,"acceleration"});  
#line 428
event_register((Event){0,1,projection,{projection_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",428,"projection"});  
#line 443
event_register((Event){0,1,end_timestep,{end_timestep_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/centered.h",443,"end_timestep"});  
#line 140 "/home/esteban/ProgramFile/basilisk/src/vof.h"
event_register((Event){0,1,stability_0,{stability_0_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/vof.h",140,"stability"});  
#line 380
event_register((Event){0,1,vof_0,{vof_0_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/vof.h",380,"vof"});  
#line 50 "/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h"
event_register((Event){0,1,tracer_advection_0,{tracer_advection_0_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h",50,"tracer_advection"});  
#line 83
event_register((Event){0,1,properties_0,{properties_0_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/two-phase-generic.h",83,"properties"});  
#line 45 "/home/esteban/ProgramFile/basilisk/src/iforce.h"
event_register((Event){0,1,acceleration_0,{acceleration_0_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/iforce.h",45,"acceleration"});  
#line 36 "/home/esteban/ProgramFile/basilisk/src/tension.h"
event_register((Event){0,1,stability_1,{stability_1_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/tension.h",36,"stability"});  
#line 72
event_register((Event){0,1,acceleration_1,{acceleration_1_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/tension.h",72,"acceleration"});  
#line 36 "/home/esteban/ProgramFile/basilisk/src/reduced.h"
event_register((Event){0,1,acceleration_2,{acceleration_2_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/reduced.h",36,"acceleration"});  
#line 72 "/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h"
event_register((Event){0,1,stability_2,{stability_2_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",72,"stability"});  
#line 117
event_register((Event){0,1,vof_1,{vof_1_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",117,"vof"});  
#line 192
event_register((Event){0,1,tracer_advection_1,{tracer_advection_1_expr0},((int *)0),((double *)0),"/home/esteban/ProgramFile/basilisk/src/navier-stokes/conserving.h",192,"tracer_advection"});
  
#line 24 "ast/init_solver.h"
}
  set_fpe();
}
