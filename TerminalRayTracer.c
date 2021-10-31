/**
 * Simple color ray tracer in the terminal
 * ray trace with the following steps
 * 1. create a scene
 * 2. create a camera
 * 3. create a ray for each pixel
 * 4. trace each ray
 * 5. display the result
 * 6. repeat
 */

//Reference frame directions
// x - right
// y - up
// z - backwards

//cameras look out from the -z direction

//lengths are in meters

//TODO->
// make pixels use colored ascii characters based on brightness
// set up the math for tracing rays around the scene
// make the scene more complex (e.g. bouncing balls)
// add proper lighting (lambertian and phong)
// add shadows? tbd on if this isn't handled already by the ray tracer algorithm

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

//math constants
#define PI 3.14159265358979323846
#define EPSILON 0.000001

//constants for the screen
#define SCREEN_WIDTH 480
#define SCREEN_HEIGHT 280

#define FRAME_RATE 60                                              //frames per second
#define FRAME_DURATION 1.0 / FRAME_RATE                            //time between frames
#define FRAME_DURATION_NS (long long)(FRAME_DURATION * 1000000000) //nanoseconds between frames

#define BOUNCE_LIMIT 10 // number of times a ray can bounce before it is considered a shadow

//TODO->come up with better method for distributing the rays uniformly within the pixel
//shoot multiple rays per pixel to handle anti-aliasing
#define RAYS_PER_PIXEL 10

//enum for the types of objects possible to hit in the scene
typedef enum
{
    NONE,
    SPHERE,
    GROUND,
    //TODO-> add more objects
} ObjectType;

//struct for a point
typedef struct
{
    double x;
    double y;
    double z;
} Point;

//struct for a vector
typedef struct
{
    double x;
    double y;
    double z;
} Vector;

Vector SKY_COLOR = {0.5372549019607842924, 0.7686274509803922017, 0.9568627450980392579};
// Vector SKY_COLOR = {0.5, 0.5, 0.85};
Vector BACKGROUND_COLOR = {0.0, 0.0, 0.0};
Vector GROUND_EVEN_COLOR = {1.0, 1.0, 1.0};
Vector GROUND_ODD_COLOR = {1.0, 0.0, 0.0};

//struct for a basis, i.e. x y and z axis of a reference frame
typedef struct
{
    Vector x;
    Vector y;
    Vector z;
} Basis;

//struct for a reference frame
typedef struct
{
    Basis basis;
    Point origin;
} Frame;

//struct for a ray
typedef struct
{
    Point origin;
    Vector direction;
} Ray;

//struct for a material (simple color only for now)
typedef struct
{
    Vector color;
    double reflectivity;
    double specularity;
} Material;

//struct for representing colors as a triple of bytes
typedef struct
{
    unsigned char r;
    unsigned char g;
    unsigned char b;
} Color;

//struct for representing a cubemap skybox. each face is a ppm texture mapped to a 1D list of colors
typedef struct
{
    Color *pX;
    Color *nX;
    Color *pY;
    Color *nY;
    Color *pZ;
    Color *nZ;
    int dim;
} Skybox;

//struct for a directional light
typedef struct
{
    Vector direction;
    Vector color;
} DirectionalLight;

//struct for a point light
typedef struct
{
    Point position;
    Vector color;
    double intensity;
} PointLight;

//struct for a sphere
typedef struct
{
    Point center;
    double radius;
    Material material;
} Sphere;

//struct for a plane
typedef struct
{
    Point point;
    Vector normal;
    Material even_material;
    Material odd_material;
} Plane;

//struct for a camera. camera is represented by a frame looking at a screen
typedef struct
{
    Frame frame;
    double screen_distance;
    double screen_width;
    double screen_height;
} Camera;

//struct for a screen consisting of a 2D array of pixels represented by vectors for color
//pixels are stored in a 1D array indexed by row*width + column
typedef struct
{
    Vector *pixels;
    int width;
    int height;
} Screen;

//struct for a scene
typedef struct
{
    Sphere *spheres;
    int num_spheres;
    Plane ground;
    DirectionalLight *directional_lights;
    int num_directional_lights;
    PointLight *point_lights;
    int num_point_lights;
    Camera camera;
    Material sky; //TODO->consider converting this to a texture map
} Scene;

//function to generate a random number between 0 and 1
double random_number()
{
    return (double)rand() / (double)RAND_MAX;
}

//function to generate a random number between min and max
double random_number_range(double min, double max)
{
    return min + random_number() * (max - min);
}

//function for generating a triangle wave pattern
//input is t.
//output: t=0->0, pi/2->1, pi->0, 3pi/2->1, 2pi->0
double triangle_wave(double t)
{
    return (fmod(t, 2 * PI) < PI) ? (fmod(t, 2 * PI) / PI) : (2 - (fmod(t, 2 * PI) / PI));
}

//function to initialize a sphere passed in with random values
void init_random_sphere(Sphere *sphere)
{
    sphere->center.x = random_number_range(-1.0, 1.0);
    sphere->center.y = random_number_range(-1.0, 1.0);
    sphere->center.z = random_number_range(-1.0, 1.0);
    sphere->radius = random_number_range(0.1, 0.5);
    sphere->material.color.x = random_number_range(0.0, 1.0);
    sphere->material.color.y = random_number_range(0.0, 1.0);
    sphere->material.color.z = random_number_range(0.0, 1.0);
}

//function to set a vector with the given values
void set_vector(Vector *vector, double x, double y, double z)
{
    vector->x = x;
    vector->y = y;
    vector->z = z;
}

//function to copy one vector into another vector
void copy_vector(Vector *dest, Vector *src)
{
    dest->x = src->x;
    dest->y = src->y;
    dest->z = src->z;
}

//function to set a point with the given values
void set_point(Point *point, double x, double y, double z)
{
    point->x = x;
    point->y = y;
    point->z = z;
}

//function to copy one point into another point
void copy_point(Point *dest, Point *src)
{
    dest->x = src->x;
    dest->y = src->y;
    dest->z = src->z;
}

//function to set a ray with the given values
void set_ray(Ray *ray, Point origin, Vector direction)
{
    ray->origin = origin;
    ray->direction = direction;
}

//function to initialize a frame with default orientation and origin
void init_frame(Frame *frame)
{
    set_vector(&(frame->basis.x), 1.0, 0.0, 0.0);
    set_vector(&(frame->basis.y), 0.0, 1.0, 0.0);
    set_vector(&(frame->basis.z), 0.0, 0.0, 1.0);
    set_point(&(frame->origin), 0.0, 0.0, 0.0);
}

//initialize a camera with default values
void init_camera(Camera *camera)
{
    init_frame(&(camera->frame));
    camera->screen_distance = 1.0;
    camera->screen_width = (double)SCREEN_WIDTH / (double)SCREEN_HEIGHT;
    camera->screen_height = 1.0;
}

//function to open a ppm file and read the data into a 1D array of colors
//takes in a Color** which is allocated with the 1D array of colors that, so the caller must free it
void read_ppm(char *filename, Color **colors_ptr, int *width, int *height)
{
    FILE *fp;

    //result return values for making fscanf and fgets not issue warnings because we ignore their returns
    char *res1;
    int res2;

    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("Error opening file %s\n", filename);
        exit(1);
    }

    //read the header
    char magic_number[3];
    res1 = fgets(magic_number, 3, fp);
    if (strncmp(magic_number, "P6", 2) != 0)
    {
        printf("Error: file is not ppm\n");
        exit(1);
    }
    fgetc(fp); //skip space/newline after the magic number

    //skip comments. comments start with a #, and end with a newline
    while (fgetc(fp) == '#')
        while (fgetc(fp) != '\n')
            ;
    fseek(fp, -1, SEEK_CUR); //back up file pointer by a single character

    //read the width and height of the image
    res2 = fscanf(fp, "%d %d", width, height);
    fgetc(fp); //skip space/newline after the width and height

    //read the maximum color value
    int max_color_value;
    res2 = fscanf(fp, "%d", &max_color_value);
    fgetc(fp); //skip space/newline after the max color value

    //assert that the max color value is 255, as we only support 8-bit colors
    if (max_color_value != 255)
    {
        printf("Error: max color value is not 255\n");
        exit(1);
    }

    //allocate memory for the colors
    *colors_ptr = (Color *)malloc(sizeof(Color) * (*width) * (*height));
    Color *colors = *colors_ptr;

    //return an error if the allocation failed
    if (colors == NULL)
    {
        printf("Error allocating memory for colors\n");
        exit(1);
    }

    //read the colors. each color is represented by 3 bytes, r, g, and b
    for (int i = 0; i < (*width) * (*height); i++)
    {
        colors[i].r = fgetc(fp);
        colors[i].g = fgetc(fp);
        colors[i].b = fgetc(fp);
    }
}

//function for loading the skybox into a Skybox struct
//skybox is selected with the name of the folder with the ppm textures in the skybox folder
//the skybox folder should contain the following files:
//  -X.ppm +X.ppm -Y.ppm +Y.ppm -Z.ppm +Z.ppm
//where X, Y, and Z are the names of the faces of the skybox
//the ppm files should be of the same size and should have matching widths and heights
void load_skybox(Skybox *skybox, char *skybox_name)
{
    char *file_names[] = {"-X.ppm", "+X.ppm", "-Y.ppm", "+Y.ppm", "-Z.ppm", "+Z.ppm"};

    // char ppm_file_name[100];
    // int ppm_width, ppm_height;
    // int ppm_size;
    // int ppm_index;
    // int ppm_row, ppm_column;
    // int ppm_color_index;
    // Color ppm_color;
    // FILE *ppm_file;

    // //load the skybox
    // sprintf(ppm_file_name, "skybox/%s/X.ppm", skybox_name);
    // ppm_file = fopen(ppm_file_name, "rb");
    // fscanf(ppm_file, "P6\n%d %d\n255\n", &ppm_width, &ppm_height);
    // ppm_size = ppm_width * ppm_height;
    // skybox->pX = (Color *)malloc(sizeof(Color) * ppm_size);
    // for (ppm_row = 0; ppm_row < ppm_height; ppm_row++)
    // {
    //     for (ppm_column = 0; ppm_column < ppm_width; ppm_column++)
    //     {
    //         ppm_color_index = ppm_row * ppm_width + ppm_column;
    //         ppm_color.r = fgetc(ppm_file);
    //         ppm_color.g = fgetc(ppm_file);
    //         ppm_color.b = fgetc(ppm_file);
    //         skybox->pX[ppm_color_index] = ppm_color;
    //     }
    // }
    // fclose(ppm_file);

    // sprintf(ppm_file_name, "skybox/%s/-X.ppm", skybox_name);
    // ppm_file = fopen(ppm_file_name, "rb");
    // fscanf(ppm_file, "P6\n%d %d\n255\n", &ppm_width, &ppm_height);
    // ppm_size = ppm_width * ppm_height;
    // skybox->mX = (Color *)malloc(sizeof(Color) * ppm_size);
    // for (ppm_row = 0; ppm_row < ppm_height; ppm_row++)
    // {
    //     for (ppm_column = 0; ppm_column < ppm_width; ppm_column++)
    //     {
    //         ppm_color_index = ppm_row * ppm_width + ppm_column;
    //         ppm_color.r = fgetc(ppm_file);
    //         ppm_color.g = fgetc(ppm_file);
    //         ppm_color.b = fgetc(ppm_file);
    //         skybox->mX[ppm_color_index] = ppm_color;
    //     }
    // }
    // fclose(ppm_file);

    // sprintf(ppm_file_name, "skybox/%s/+X.ppm", skybox_name);
    // ppm_file = fopen(ppm_file_name, "rb");
    // fscanf(ppm_file, "P6\n%d %d\n255\n", &ppm_width, &ppm_height);
    // ppm_size = ppm_width * ppm_height;
    // skybox->pX = (Color *)malloc(sizeof(Color) * ppm_size);
    // for (ppm_row = 0; ppm_row < ppm_height; ppm_row++)
    // {
    //     for (ppm_column = 0; ppm_column < ppm_width; ppm_column++)
    //     {
    //         ppm_color_index = ppm_row * ppm_width + ppm_column;
    //         ppm_color.r = fgetc(ppm_file);
    //         ppm_color.g = fgetc(ppm_file);
    //         ppm_color.b = fgetc(ppm_file);
    //         skybox->pX[ppm_color_index] = ppm_color;
    //     }
    // }
}

//function to normalize a vector
void normalize_vector(Vector *vector)
{
    double length = sqrt(vector->x * vector->x + vector->y * vector->y + vector->z * vector->z);

    // only normalize if length isn't basically 0
    if (length > 0.0001)
    {
        vector->x /= length;
        vector->y /= length;
        vector->z /= length;
    }
}

//compute the dot product of two vectors
double dot_product(Vector *vector1, Vector *vector2)
{
    return vector1->x * vector2->x + vector1->y * vector2->y + vector1->z * vector2->z;
}

//multiply a vector by a scalar
void scale_vector(Vector *vector, double scalar)
{
    vector->x *= scalar;
    vector->y *= scalar;
    vector->z *= scalar;
}

//return a scaled vector without modifying the original
Vector scale_vector_copy(Vector *vector, double scalar)
{
    Vector result = *vector;
    scale_vector(&result, scalar);
    return result;
}

//add two vectors and store the result in vector 1
void add_vectors(Vector *vector1, Vector *vector2)
{
    vector1->x += vector2->x;
    vector1->y += vector2->y;
    vector1->z += vector2->z;
}

//add two vectors without modifying either, returning the result
Vector add_vectors_copy(Vector *vector1, Vector *vector2)
{
    Vector result = *vector1;
    add_vectors(&result, vector2);
    return result;
}

//subtract vector 2 from vector 1 and store the result in vector 1
void subtract_vectors(Vector *vector1, Vector *vector2)
{
    vector1->x -= vector2->x;
    vector1->y -= vector2->y;
    vector1->z -= vector2->z;
}

//subtract vector 2 from vector 1 without modifying either, returning the result
Vector subtract_vectors_copy(Vector *vector1, Vector *vector2)
{
    Vector result = *vector1;
    subtract_vectors(&result, vector2);
    return result;
}

//function to multiply elements of two vectors point-wise
void multiply_vectors(Vector *vector1, Vector *vector2)
{
    vector1->x *= vector2->x;
    vector1->y *= vector2->y;
    vector1->z *= vector2->z;
}

//function to clamp a number within the given range
double clamp(double value, double min, double max)
{
    if (value < min)
        return min;
    if (value > max)
        return max;
    return value;
}

//function to clamp parameters of a vector to the given range
void clamp_vector(Vector *vector, double min, double max)
{
    vector->x = clamp(vector->x, min, max);
    vector->y = clamp(vector->y, min, max);
    vector->z = clamp(vector->z, min, max);
}

//function to multiply elements of two vectors point-wise without modifying either
Vector multiply_vectors_copy(Vector *vector1, Vector *vector2)
{
    Vector result = *vector1;
    multiply_vectors(&result, vector2);
    return result;
}

// //function to compute the distance between two points
// double point_distance(Point *point1, Point *point2)
// {
//     double dx = point1->x - point2->x;
//     double dy = point1->y - point2->y;
//     double dz = point1->z - point2->z;
//     return sqrt(dx * dx + dy * dy + dz * dz);
// }

//function to rotate a basis by another basis (rotation matrix) and ensure orthonomal
void rotate_basis(Basis *basis, Basis *rotation)
{
    Basis result;
    result.x = (Vector){.x = basis->x.x * rotation->x.x + basis->x.y * rotation->x.y + basis->x.z * rotation->x.z,
                        .y = basis->x.x * rotation->y.x + basis->x.y * rotation->y.y + basis->x.z * rotation->y.z,
                        .z = basis->x.x * rotation->z.x + basis->x.y * rotation->z.y + basis->x.z * rotation->z.z};
    result.y = (Vector){.x = basis->y.x * rotation->x.x + basis->y.y * rotation->x.y + basis->y.z * rotation->x.z,
                        .y = basis->y.x * rotation->y.x + basis->y.y * rotation->y.y + basis->y.z * rotation->y.z,
                        .z = basis->y.x * rotation->z.x + basis->y.y * rotation->z.y + basis->y.z * rotation->z.z};
    result.z = (Vector){.x = basis->z.x * rotation->x.x + basis->z.y * rotation->x.y + basis->z.z * rotation->x.z,
                        .y = basis->z.x * rotation->y.x + basis->z.y * rotation->y.y + basis->z.z * rotation->y.z,
                        .z = basis->z.x * rotation->z.x + basis->z.y * rotation->z.y + basis->z.z * rotation->z.z};
    basis->x = result.x;
    basis->y = result.y;
    basis->z = result.z;
}

//function to rotate a basis by the given angle around the x axis
void rotate_basis_x(Basis *basis, double angle)
{
    Basis rotation;
    set_vector(&(rotation.x), 1.0, 0.0, 0.0);
    set_vector(&(rotation.y), 0.0, cos(angle), -sin(angle));
    set_vector(&(rotation.z), 0.0, sin(angle), cos(angle));
    rotate_basis(basis, &rotation);
}

//function to rotate a basis by the given angle around the y axis
void rotate_basis_y(Basis *basis, double angle)
{
    Basis rotation;
    set_vector(&(rotation.x), cos(angle), 0.0, sin(angle));
    set_vector(&(rotation.y), 0.0, 1.0, 0.0);
    set_vector(&(rotation.z), -sin(angle), 0.0, cos(angle));
    rotate_basis(basis, &rotation);
}

//function to rotate a basis by the given angle around the z axis
void rotate_basis_z(Basis *basis, double angle)
{
    Basis rotation;
    set_vector(&(rotation.x), cos(angle), -sin(angle), 0.0);
    set_vector(&(rotation.y), sin(angle), cos(angle), 0.0);
    set_vector(&(rotation.z), 0.0, 0.0, 1.0);
    rotate_basis(basis, &rotation);
}

//function to transform a reference frame by the position and rotation of another frame
//performs the 4x4 matrix transform treating the basis and position as a single 4x4 homogeneous matrix
void transform_frame(Frame *frame, Frame *transform)
{
    Basis result_basis = {
        .x = {.x = frame->basis.x.x * transform->basis.x.x + frame->basis.x.y * transform->basis.y.x + frame->basis.x.z * transform->basis.z.x,
              .y = frame->basis.x.x * transform->basis.x.y + frame->basis.x.y * transform->basis.y.y + frame->basis.x.z * transform->basis.z.y,
              .z = frame->basis.x.x * transform->basis.x.z + frame->basis.x.y * transform->basis.y.z + frame->basis.x.z * transform->basis.z.z},
        .y = {.x = frame->basis.y.x * transform->basis.x.x + frame->basis.y.y * transform->basis.y.x + frame->basis.y.z * transform->basis.z.x,
              .y = frame->basis.y.x * transform->basis.x.y + frame->basis.y.y * transform->basis.y.y + frame->basis.y.z * transform->basis.z.y,
              .z = frame->basis.y.x * transform->basis.x.z + frame->basis.y.y * transform->basis.y.z + frame->basis.y.z * transform->basis.z.z},
        .z = {.x = frame->basis.z.x * transform->basis.x.x + frame->basis.z.y * transform->basis.y.x + frame->basis.z.z * transform->basis.z.x,
              .y = frame->basis.z.x * transform->basis.x.y + frame->basis.z.y * transform->basis.y.y + frame->basis.z.z * transform->basis.z.y,
              .z = frame->basis.z.x * transform->basis.x.z + frame->basis.z.y * transform->basis.y.z + frame->basis.z.z * transform->basis.z.z}};
    Point result_position = {.x = frame->origin.x * transform->basis.x.x + frame->origin.y * transform->basis.y.x + frame->origin.z * transform->basis.z.x + transform->origin.x,
                             .y = frame->origin.x * transform->basis.x.y + frame->origin.y * transform->basis.y.y + frame->origin.z * transform->basis.z.y + transform->origin.y,
                             .z = frame->origin.x * transform->basis.x.z + frame->origin.y * transform->basis.y.z + frame->origin.z * transform->basis.z.z + transform->origin.z};
    frame->basis = result_basis;
    frame->origin = result_position;
}

//function to reflect a vector about a normal vector
void reflect_vector(Vector *vector, Vector *normal)
{
    double dot = dot_product(vector, normal);
    vector->x = vector->x - 2.0 * dot * normal->x;
    vector->y = vector->y - 2.0 * dot * normal->y;
    vector->z = vector->z - 2.0 * dot * normal->z;
}

//function to compute the location of a ray intersection with a sphere
//returns 1 if the ray intersects, 0 otherwise
//intersection point is stored in the provided reference to an output point
int ray_intersects_sphere(Ray *ray, Sphere *sphere, Point *intersection_point)
{
    Vector ray_origin_sphere_vector;
    set_vector(&ray_origin_sphere_vector,
               ray->origin.x - sphere->center.x,
               ray->origin.y - sphere->center.y,
               ray->origin.z - sphere->center.z);

    double a = dot_product(&(ray->direction), &(ray->direction));
    double b = 2.0 * dot_product(&ray_origin_sphere_vector, &(ray->direction));
    double c = dot_product(&ray_origin_sphere_vector, &ray_origin_sphere_vector) - sphere->radius * sphere->radius;

    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0)
    {
        return 0;
    }
    else
    {
        double t0 = (-b - sqrt(discriminant)) / (2.0 * a);
        // double t1 = (-b + sqrt(discriminant)) / (2.0 * a);
        if (t0 > 0.0)
        {
            set_vector((Vector *)intersection_point,
                       ray->origin.x + t0 * ray->direction.x,
                       ray->origin.y + t0 * ray->direction.y,
                       ray->origin.z + t0 * ray->direction.z);
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

//function to compute the location of a ray intersection with a plane
//returns 1 if the ray intersects, 0 otherwise
//intersection point is stored in the provided reference to an output point
int ray_intersects_plane(Ray *ray, Plane *plane, Point *intersection_point)
{
    double denom = dot_product(&(ray->direction), &(plane->normal));
    if (fabs(denom) > 0.00001)
    {
        // Vector negative_origin = scale_vector_copy(&ray->origin, -1.0);
        Vector ray_to_play_point = subtract_vectors_copy((Vector *)&plane->point, (Vector *)&ray->origin);
        double t = dot_product(&ray_to_play_point, &(plane->normal)) / denom;
        if (t > 0.00001)
        {
            set_vector((Vector *)intersection_point,
                       ray->origin.x + t * ray->direction.x,
                       ray->origin.y + t * ray->direction.y,
                       ray->origin.z + t * ray->direction.z);
            return 1;
        }
    }
    return 0;
}

//break out function for computing the closest intersection point of an object in the scene
//returns closest intersection point, normal of the intersection, and material of intersection point
ObjectType trace_ray(Scene *scene, Ray *ray, Point *intersection, Vector *normal, Material *material)
{
    //find the closest intersection
    double closest_distance = INFINITY;
    ObjectType closest_object = NONE;
    int closest_index = -1;
    Point point;                //point to use when computing the intersection point
    Point closest_intersection; //closest intersection point
    Vector closest_normal;      //normal of the closest intersection point
    Material closest_material;  //material of the closest intersection point

    //check each of the spheres if they are the closest intersection
    for (int i = 0; i < scene->num_spheres; i++)
    {
        if (ray_intersects_sphere(ray, &(scene->spheres[i]), &point))
        {
            //efficiently compute distance from ray origin to intersection point
            Vector ray_origin_sphere_vector;
            set_vector(&ray_origin_sphere_vector,
                       ray->origin.x - point.x,
                       ray->origin.y - point.y,
                       ray->origin.z - point.z);
            double square_distance = dot_product(&ray_origin_sphere_vector, &ray_origin_sphere_vector);
            if (square_distance < closest_distance)
            {
                closest_object = SPHERE;
                closest_distance = square_distance;
                closest_index = i;

                //assign the output values
                closest_intersection = point;
                closest_normal = subtract_vectors_copy((Vector *)&point, (Vector *)&(scene->spheres[i].center));
                closest_material = scene->spheres[i].material;
            }
        }
    }

    //check the ground if it is the closest intersection
    if (ray_intersects_plane(ray, &(scene->ground), &point))
    {
        //efficiently compute distance from ray origin to intersection point
        Vector ray_origin_plane_vector;
        set_vector(&ray_origin_plane_vector,
                   ray->origin.x - point.x,
                   ray->origin.y - point.y,
                   ray->origin.z - point.z);
        double square_distance = dot_product(&ray_origin_plane_vector, &ray_origin_plane_vector);
        if (square_distance < closest_distance)
        {
            closest_object = GROUND;
            closest_distance = square_distance;

            //assign the output values if they aren't null
            closest_intersection = point;
            closest_normal = scene->ground.normal;

            //ground material is a checker pattern based on if the intersection point's coordinates are even/odd
            int odd = (int)(floor(point.x) + floor(point.z)) & 1;
            closest_material = odd ? scene->ground.odd_material : scene->ground.even_material;
        }
    }

    //TODO->check other object types

    //if no object was intersected, set the return values to defaults, and the sky color
    if (closest_object == NONE)
    {
        closest_intersection = ray->origin;
        closest_normal = ray->direction;
        closest_material = scene->sky;
    }
    else
    {
        //push the intersection point back a little bit to avoid self-intersection
        Vector to_surface = subtract_vectors_copy((Vector *)&ray->origin, (Vector *)&closest_intersection);
        normalize_vector(&to_surface);
        scale_vector(&to_surface, EPSILON);
        add_vectors((Vector *)&closest_intersection, &to_surface);
    }

    //ensure the normal is normalized
    normalize_vector(&closest_normal);

    //return the closest intersection point, normal, and material if they are not null
    if (intersection != NULL)
        *intersection = closest_intersection;
    if (normal != NULL)
        *normal = closest_normal;
    if (material != NULL)
        *material = closest_material;

    return closest_object;
}

//function to apply lighting to a point in the scene. resulting color is stored in color which contains the initial unlit color of the object
//specular highlights are computed using the Blinn-Phong model
//view points from the object to the source of the ray for this point
void apply_lighting(Scene *scene, Point *intersection, Vector *view, Vector *normal, Material *material)
{
    //vector to accumulate the color of point based on the lighting
    Vector output_color = {0.0, 0.0, 0.0};

    //compute the contribution of the directional lights
    for (int i = 0; i < scene->num_directional_lights; i++)
    {
        //create a ray for this light source to see if it is blocked by any objects in the scene
        Vector light_direction = scale_vector_copy(&(scene->directional_lights[i].direction), -1.0);
        normalize_vector(&light_direction);

        //check if the light is blocked by any objects in the scene
        ObjectType blocking_object = trace_ray(scene, &(Ray){.origin = *intersection, .direction = light_direction}, NULL, NULL, NULL);
        if (blocking_object == NONE)
        {
            //compute the diffuse and specular contributions
            Vector diffuse_contribution = scale_vector_copy(&scene->directional_lights[i].color, fmin(dot_product(normal, &light_direction), 1.0));

            // //compute the blinn-phong contribution
            // Vector half = add_vectors_copy(&light_direction, view);
            // normalize_vector(&half);
            // Vector specular_contribution = scale_vector_copy(&scene->directional_lights[i].color, pow(clamp(dot_product(normal, &half), 0.0, 1.0), material->specularity));

            //add the contributions to the color
            multiply_vectors(&diffuse_contribution, &material->color);
            add_vectors(&output_color, &diffuse_contribution);
            // add_vectors(&output_color, &specular_contribution);
        }
    }

    //compute the contribution of the point lights
    for (int i = 0; i < scene->num_point_lights; i++)
    {
        //create a ray for this light source to see if it is blocked by any objects in the scene
        Vector light_direction = subtract_vectors_copy((Vector *)&(scene->point_lights[i].position), (Vector *)intersection);
        double light_distance_squared = dot_product(&light_direction, &light_direction);
        double light_intensity = clamp(scene->point_lights[i].intensity / light_distance_squared, 0.0, 1.0);
        // light_distance = sqrt(light_distance);
        normalize_vector(&light_direction);

        //check if the light is blocked by any objects in the scene
        Point point;
        ObjectType blocking_object = trace_ray(scene, &(Ray){.origin = *intersection, .direction = light_direction}, &point, NULL, NULL);
        //distance to intersection
        Vector to_intersection = subtract_vectors_copy((Vector *)&point, (Vector *)intersection);
        double intersection_distance_squared = dot_product(&to_intersection, &to_intersection);

        if (blocking_object == NONE || light_distance_squared < intersection_distance_squared)
        {
            //compute the diffuse and specular contributions
            Vector diffuse_contribution = scale_vector_copy(&scene->point_lights[i].color, light_intensity * fmin(dot_product(normal, &light_direction), 1.0));

            // //compute the blinn-phong contribution
            // Vector half = add_vectors_copy(&light_direction, view);
            // normalize_vector(&half);
            // Vector specular_contribution = scale_vector_copy(&scene->point_lights[i].color, light_intensity * pow(clamp(dot_product(normal, &half), 0.0, 1.0), material->specularity));

            //add the contributions to the color
            multiply_vectors(&diffuse_contribution, &material->color);
            add_vectors(&output_color, &diffuse_contribution);
            // add_vectors(&output_color, &specular_contribution);
        }
    }

    //clamp the color so that it is not greater than 1.0
    clamp_vector(&output_color, 0.0, 1.0);

    material->color = output_color;
}

//function to project the scene onto a screen
void project_scene(Scene *scene, Screen *screen)
{
    //get the camera frame y and z vectors
    Vector camera_frame_x = scene->camera.frame.basis.x;
    Vector camera_frame_y = scene->camera.frame.basis.y;

    //iterate over each pixel in the screen
    for (int row = 0; row < screen->height; row++)
    {
        for (int column = 0; column < screen->width; column++)
        {
            Vector average_pixel_color = {.x = 0.0, .y = 0.0, .z = 0.0};
            for (int ray_num = 0; ray_num < RAYS_PER_PIXEL; ray_num++)
            {
                //compute the size of a pixel on the screen
                double pixel_width = scene->camera.screen_width / screen->width;
                double pixel_height = scene->camera.screen_height / screen->height;

                //compute the ray for the pixel

                //point on the screen in 2D space
                double screen_x = (((double)column / (double)screen->width) * scene->camera.screen_width - scene->camera.screen_width / 2.0);
                double screen_y = -(((double)row / (double)screen->height) * scene->camera.screen_height - scene->camera.screen_height / 2.0);
                double screen_z = -scene->camera.screen_distance;

                //offset the x and y coordinates by the sub pixel offset for antialiasing
                screen_x += triangle_wave(2 * PI * ray_num / RAYS_PER_PIXEL) / 2 * pixel_width;
                screen_y += triangle_wave(PI * ray_num / RAYS_PER_PIXEL) / 2 * pixel_height;

                //convert point to 3D space relative to the scene origin
                Vector world_x = scale_vector_copy(&scene->camera.frame.basis.x, screen_x);
                Vector world_y = scale_vector_copy(&scene->camera.frame.basis.y, screen_y);
                Vector world_z = scale_vector_copy(&scene->camera.frame.basis.z, screen_z);
                Vector screen_point = (Vector){0.0, 0.0, 0.0};
                add_vectors(&screen_point, &world_x);
                add_vectors(&screen_point, &world_y);
                add_vectors(&screen_point, &world_z);

                //make the point point from the camera origin to the screen point
                subtract_vectors(&screen_point, (Vector *)&scene->camera.frame.origin);

                //normalize the ray direction
                normalize_vector(&screen_point);

                //create the ray from the screen_point
                Ray ray = (Ray){scene->camera.frame.origin, screen_point};
                Vector pixel_color = {.x = 0.0, .y = 0.0, .z = 0.0};
                int bounces = 0;
                double color_contribution = 1.0;       //amount of color to the pixel for the current bounce. updated by material reflectivity
                double color_contribution_total = 0.0; //total amount of color contributed by all bounces
                int still_going = 1;                   //false if the ray has stoped bouncing off of things

                while (still_going && bounces < BOUNCE_LIMIT && color_contribution > 0.00001)
                {
                    // shoot the ray and find the closest intersection
                    Point intersection;
                    Vector normal;
                    Material material;
                    ObjectType closest_object = trace_ray(scene, &ray, &intersection, &normal, &material);

                    if (closest_object != NONE)
                    {
                        //determine the apparent color of the intersection point based on the lighting in the scene
                        Vector view = scale_vector_copy(&ray.direction, -1.0);
                        apply_lighting(scene, &intersection, &view, &normal, &material);
                    }

                    //accumulate the total color contribution, and update the color to merge into the pixel according to the contribution
                    color_contribution_total += color_contribution;
                    scale_vector(&material.color, color_contribution);

                    //if hit an object, adjust the new color contribution based on the material reflectivity and the amount of bounces
                    //else stop the ray
                    if (closest_object != NONE)
                    {
                        color_contribution *= material.reflectivity;
                        bounces++;
                    }
                    else //hit the sky
                    {
                        color_contribution = 0.0;
                        still_going = 0;
                    }

                    //add the color contribution from this bounce
                    add_vectors(&pixel_color, &material.color);

                    //compute the reflection ray
                    reflect_vector(&ray.direction, &normal);
                    normalize_vector(&ray.direction);
                    copy_point(&ray.origin, &intersection);
                }

                //set the color of this pixel in the screen

                scale_vector(&pixel_color, 1.0 / color_contribution_total);

                add_vectors(&average_pixel_color, &pixel_color);
            }
            scale_vector(&average_pixel_color, 1.0 / RAYS_PER_PIXEL);
            screen->pixels[row * screen->width + column] = average_pixel_color;
        }
    }
}

//function to set all pixels in the screen to a given color
void set_screen_color(Screen *screen, Vector color)
{
    for (int row = 0; row < screen->height; row++)
    {
        for (int column = 0; column < screen->width; column++)
        {
            screen->pixels[row * screen->width + column] = color;
        }
    }
}

//function to draw a screen to the terminal
void draw_screen(Screen *screen)
{
    //move the cursor back to the top of the terminal
    printf("\033[0;0H");

    //print out each pixel in the screen
    for (int i = 0; i < screen->height; i++)
    {
        for (int j = 0; j < screen->width; j++)
        {
            Vector pixel = screen->pixels[i * screen->width + j];
            printf("\033[48;2;%d;%d;%dm  \033[0m", (int)(pixel.x * 255), (int)(pixel.y * 255), (int)(pixel.z * 255));
        }
        printf("\n");
    }
}

//more efficient function to draw to the screen using a string buffer and single print statement
char reset_str[] = "\033[0;0H";
char pixel_str[] = "\033[48;2;000;000;000m  \033[0m";
char screenbuffer[(sizeof(reset_str) + 1) + ((sizeof(pixel_str) - 1) * SCREEN_WIDTH + 1) * SCREEN_HEIGHT + 1];

//function to initialize screenbuffer with color codes at every pixel location
void initialize_screenbuffer()
{
    char *buff_ptr = screenbuffer;

    //copy the reset string to the front of the buffer, excluding the null terminator
    memcpy(buff_ptr, reset_str, sizeof(reset_str) - 1);

    //increment the buffer pointer past the reset string
    buff_ptr += sizeof(reset_str) - 1;

    //copy the pixel string to the rest of the buffer, excluding the null terminator, and adding a newline after each row
    for (int i = 0; i < SCREEN_HEIGHT; i++)
    {
        for (int j = 0; j < SCREEN_WIDTH; j++)
        {
            memcpy(buff_ptr, pixel_str, sizeof(pixel_str) - 1);
            buff_ptr += sizeof(pixel_str) - 1;
        }
        memcpy(buff_ptr, "\n", 1);
        buff_ptr += 1;
    }

    //null terminate the buffer
    *buff_ptr = '\0';
}

//function to convert a byte to its 3 base-10 digit values
void byte_to_digits(int value, char output[3])
{
    output[0] = value / 100 + '0';
    output[1] = (value / 10) % 10 + '0';
    output[2] = value % 10 + '0';
}

//function to draw screen by setting the values in the screenbuffer and printing it in one print statement
void buffered_draw_screen(Screen *screen)
{
    char *buff_ptr = screenbuffer;

    //move the ptr past the reset string to the first row of pixels
    buff_ptr += sizeof(reset_str) - 1;

    //set the value of each pixel color in the buffer
    for (int i = 0; i < screen->height; i++)
    {
        for (int j = 0; j < screen->width; j++)
        {
            buff_ptr += 7; //advance the pointer past the start of the color code sequence
            Vector pixel = screen->pixels[i * screen->width + j];
            char val[3]; //holds the value of each color component from 0-255 (with leading 0s so it's always 3 characters long)
            byte_to_digits((int)(pixel.x * 255), val);
            memcpy(buff_ptr, val, 3);
            buff_ptr += 4; //include the semicolon after the color value
            byte_to_digits((int)(pixel.y * 255), val);
            memcpy(buff_ptr, val, 3);
            buff_ptr += 4; //include the semicolen after the color value
            byte_to_digits((int)(pixel.z * 255), val);
            memcpy(buff_ptr, val, 3);
            buff_ptr += 10; //inclide the "m  \033[0m" characters after the color value
        }
        buff_ptr += 1; //advance past the newline
    }

    //print the screenbuffer to the terminal as efficiently as possible
    fwrite(screenbuffer, sizeof(char), sizeof(screenbuffer), stdout);
}

//function to print a ppm file to the terminal as an image for testing
void print_ppm(Color *colors, int width, int height)
{
    //convert the colors to a Screen object
    Screen screen;
    screen.width = width;
    screen.height = height;
    screen.pixels = (Vector *)malloc(sizeof(Vector) * width * height);
    for (int i = 0; i < width * height; i++)
    {
        screen.pixels[i].x = colors[i].r / 255.0;
        screen.pixels[i].y = colors[i].g / 255.0;
        screen.pixels[i].z = colors[i].b / 255.0;
    }

    //print the image to the terminal
    draw_screen(&screen);

    //free the memory allocated for the pixels
    free(screen.pixels);
}

//function to capture user input arrow keys for camera movement using getch()
// void get_camera_movement(Scene *scene)
// {
//     int key = getch();
//     switch (key)
//     {
//     case KEY_UP:
//         scene->camera.frame.origin.y += CAMERA_MOVE_SPEED;
//         break;
//     case KEY_DOWN:
//         scene->camera.frame.origin.y -= CAMERA_MOVE_SPEED;
//         break;
//     case KEY_LEFT:
//         scene->camera.frame.origin.x -= CAMERA_MOVE_SPEED;
//         break;
//     case KEY_RIGHT:
//         scene->camera.frame.origin.x += CAMERA_MOVE_SPEED;
//         break;
//     }
// }

//test/create a simple scene with a single sphere in the middle and the camera looking directly at it
int main()
{
    //seed the random number generator
    srand(time(0));

    //initialize the screen buffer
    initialize_screenbuffer();

    //get the time that the program starts
    struct timespec start;
    timespec_get(&start, TIME_UTC);

    //create a list of 6 spheres 1 for each direction in 3D
    // #define NUM_SPHERES 6
    //objects in the scene
    Sphere spheres[] = {
        {.center = {.x = 0.25, .y = 0.0, .z = 0.0}, .material = {.color = {.x = 1.0, .y = 0.0, .z = 0.0}, .reflectivity = 1.0, .specularity = 100.0}, .radius = 0.125},
        {.center = {.x = 0.0, .y = 0.25, .z = 0.0}, .material = {.color = {.x = 0.0, .y = 1.0, .z = 0.0}, .reflectivity = 0.8, .specularity = 100.0}, .radius = 0.125},
        {.center = {.x = 0.0, .y = 0.0, .z = 0.25}, .material = {.color = {.x = 0.0, .y = 0.0, .z = 1.0}, .reflectivity = 0.8, .specularity = 100.0}, .radius = 0.125},
        {.center = {.x = -0.25, .y = 0.0, .z = 0.0}, .material = {.color = {.x = 0.0, .y = 1.0, .z = 1.0}, .reflectivity = 0.8, .specularity = 100.0}, .radius = 0.125},
        {.center = {.x = 0.0, .y = -0.25, .z = 0.0}, .material = {.color = {.x = 1.0, .y = 0.0, .z = 1.0}, .reflectivity = 0.8, .specularity = 100.0}, .radius = 0.125},
        {.center = {.x = 0.0, .y = 0.0, .z = -0.25}, .material = {.color = {.x = 1.0, .y = 1.0, .z = 0.0}, .reflectivity = 0.8, .specularity = 100.0}, .radius = 0.125},
    };
    const int NUM_SPHERES = sizeof(spheres) / sizeof(Sphere);

    Plane ground = {
        .normal = {.x = 0.0, .y = 1.0, .z = 0.0},
        .point = {.x = 0.0, .y = -2.0, .z = 0.0},
        .even_material = {.color = GROUND_EVEN_COLOR, .reflectivity = 0.2, .specularity = 100.0},
        .odd_material = {.color = GROUND_ODD_COLOR, .reflectivity = 0.2, .specularity = 100.0},
    };
    //TODO->other objects

    //lights in the scene
    DirectionalLight directional_lights[] = {{
        .direction = {.x = -1.0, .y = -1.0, .z = -1.0},
        .color = {.x = 1.0, .y = 1.0, .z = 1.0},
    }};
    const int NUM_DIRECTIONAL_LIGHTS = sizeof(directional_lights) / sizeof(DirectionalLight);
    PointLight point_lights[] = {
        {.position = {.x = 0.0, .y = 0.0, .z = 0.0}, .color = {.x = 1.0, .y = 1.0, .z = 1.0}, .intensity = 10.0},
        // {.position = {.x = 0.0, .y = 0.0, .z = 0.0}, .color = {.x = 1.0, .y = 1.0, .z = 1.0}, .intensity = 1.0},
        // {.position = {.x = 0.0, .y = 0.0, .z = 0.0}, .color = {.x = 1.0, .y = 1.0, .z = 1.0}, .intensity = 1.0},
    };
    const int NUM_POINT_LIGHTS = sizeof(point_lights) / sizeof(PointLight);

    //create a camera looking at the sphere from 2 meters away
    //camera position will be set inside the while loop
    Camera camera;
    init_camera(&camera);

    //create a scene with the camera and the spheres
    Scene scene = {
        .camera = camera,
        .spheres = spheres,
        .num_spheres = NUM_SPHERES,
        .ground = ground,
        .directional_lights = directional_lights,
        .num_directional_lights = NUM_DIRECTIONAL_LIGHTS,
        .point_lights = point_lights,
        .num_point_lights = NUM_POINT_LIGHTS,
        .sky = {.color = SKY_COLOR, .reflectivity = 0.0}};

    //pointer to the camera so we can update it every frame
    Camera *camera_ptr = &scene.camera;

    //create a screen of size SCREEN_WIDTH x SCREEN_HEIGHT with statically allocated pixels
    Vector pixels[SCREEN_HEIGHT * SCREEN_WIDTH];
    Screen screen = {.pixels = (Vector *)pixels, .width = SCREEN_WIDTH, .height = SCREEN_HEIGHT};

    //loop forever. set each pixel on the screen to a random color and draw the screen to the terminal
    struct timespec ts; //keep track of the system time for computing the duration of frames
    while (1)
    {
        //get the current system time, get the number of nanoseconds elapsed since the start of the program
        timespec_get(&ts, TIME_UTC);
        long long start_nanos = (ts.tv_sec - start.tv_sec) * 1000000000 + ts.tv_nsec - start.tv_nsec;

        //seconds that this frame occurred at
        double t = (double)start_nanos / 1000000000.0;

        //construct the transform of the camera. the camera should orbit the center of the scene
        Frame tf0, tf1;
        init_frame(&tf0);
        init_frame(&tf1);
        init_frame(&(scene.camera.frame));
        rotate_basis_x(&tf0.basis, 2.0 * PI * t * -0.005);
        rotate_basis_y(&tf0.basis, 2.0 * PI * t * 0.003);
        Vector root_to_camera = {.x = 0.0, .y = 0.0, .z = 1.99};
        add_vectors((Vector *)&tf1.origin, &root_to_camera);
        transform_frame(&scene.camera.frame, &tf1);
        transform_frame(&scene.camera.frame, &tf0);

        //project the scene onto the screen
        project_scene(&scene, &screen);

        //draw the screen to the terminal
        buffered_draw_screen(&screen);

        //compute the amount of time the frame computations took
        timespec_get(&ts, TIME_UTC);
        long long end_nanos = (ts.tv_sec - start.tv_sec) * 1000000000 + ts.tv_nsec - start.tv_nsec;
        long long frame_time_nanos = end_nanos - start_nanos;

        //sleep the remaining time in the frame
        if (FRAME_DURATION_NS > frame_time_nanos)
        {
            long long nanos_to_sleep = FRAME_DURATION_NS - frame_time_nanos;
            const struct timespec delay = {.tv_sec = nanos_to_sleep / 1000000000, .tv_nsec = nanos_to_sleep % 1000000000};
            nanosleep(&delay, NULL);
        }

        //compute the duration of the frame including any sleep time
        timespec_get(&ts, TIME_UTC);
        long long end_nanos_2 = (ts.tv_sec - start.tv_sec) * 1000000000 + ts.tv_nsec - start.tv_nsec;
        long long frame_duration_nanos = end_nanos_2 - start_nanos;

        //print out the frame rate at the top of the screen, and move the cursor back
        fputs("\033[0;0H", stdout);
        printf("%.02f fps\n", 1.0 / ((double)frame_duration_nanos / 1000000000.0));
        fputs("\033[0;0H", stdout);
    }
}