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

//math constants
#define PI 3.14159265358979323846

//constants for the screen
#define SCREEN_WIDTH 80
#define SCREEN_HEIGHT 80

#define FRAME_RATE 60                                              //frames per second
#define FRAME_DURATION 1.0 / FRAME_RATE                            //time between frames
#define FRAME_DURATION_NS (long long)(FRAME_DURATION * 1000000000) //nanoseconds between frames

#define BOUNCE_LIMIT 10 // number of times a ray can bounce before it is considered a shadow

//constants for the scene
#define GROUND_PLANE_HEIGHT -2.0

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

Vector SKY_COLOR = {137.0 / 255, 196.0 / 255, 244.0 / 255};
Vector BACKGROUND_COLOR = {0.0, 0.0, 0.0};

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
} Material;

//struct for a sphere
typedef struct
{
    Point center;
    double radius;
    Material material;
} Sphere;

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
    Camera camera;
    Screen screen;
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

//function to set a point with the given values
void set_point(Point *point, double x, double y, double z)
{
    point->x = x;
    point->y = y;
    point->z = z;
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
    camera->screen_width = 1.0;
    camera->screen_height = 1.0;
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

// //function to orbit the camera around the origin
// void orbit_camera(Camera *camera, double angle)
// {
//     rotate_basis_y(&(camera->frame.basis), angle);
// }

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

        // if (t1 > 0.0)
        // {
        //     set_vector((Vector *)intersection_point,
        //                ray->origin.x + ray->direction.x * t1,
        //                ray->origin.y + ray->direction.y * t1,
        //                ray->origin.z + ray->direction.z * t1);
        //     return 1;
        // }
        // else
        // {
        //     return 0;
        // }
        // if (t1 < 0.0 && t2 < 0.0)
        // {
        //     return 0;
        // }
        // else if (t1 < 0.0)
        // {
        //     set_vector((Vector *)intersection_point,
        //                ray->origin.x + ray->direction.x * t2,
        //                ray->origin.y + ray->direction.y * t2,
        //                ray->origin.z + ray->direction.z * t2);
        // }
        // else if (t2 < 0.0)
        // {
        //     set_vector((Vector *)intersection_point,
        //                ray->origin.x + ray->direction.x * t1,
        //                ray->origin.y + ray->direction.y * t1,
        //                ray->origin.z + ray->direction.z * t1);
        // }
        // else
        // {
        //     double t = fmin(t1, t2);
        //     set_vector((Vector *)intersection_point,
        //                ray->origin.x + ray->direction.x * t,
        //                ray->origin.y + ray->direction.y * t,
        //                ray->origin.z + ray->direction.z * t);
        // }
        // return 1;
    }
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
            //compute the ray for the pixel

            //point on the screen in 2D space
            double screen_x = (((double)column / (double)screen->width) * scene->camera.screen_width - scene->camera.screen_width / 2.0);
            double screen_y = -(((double)row / (double)screen->height) * scene->camera.screen_height - scene->camera.screen_height / 2.0);
            double screen_z = -scene->camera.screen_distance;

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
            int still_going = 1; //false if the ray has stoped bouncing off of things

            while (still_going && bounces < BOUNCE_LIMIT)
            {
                //find the closest intersection
                double closest_distance = INFINITY;
                int closest_sphere_index = -1;
                Point closest_intersection_point = (Point){0.0, 0.0, 0.0};
                for (int i = 0; i < scene->num_spheres; i++)
                {
                    Point intersection_point;
                    if (ray_intersects_sphere(&ray, &(scene->spheres[i]), &intersection_point))
                    {
                        //efficiently compute distance from ray origin to intersection point
                        Vector ray_origin_sphere_vector;
                        set_vector(&ray_origin_sphere_vector,
                                   ray.origin.x - intersection_point.x,  //scene->spheres[i].center.x,
                                   ray.origin.y - intersection_point.y,  //scene->spheres[i].center.y,
                                   ray.origin.z - intersection_point.z); //scene->spheres[i].center.z);
                        double square_distance = dot_product(&ray_origin_sphere_vector, &ray_origin_sphere_vector);
                        if (square_distance < closest_distance)
                        {
                            closest_distance = square_distance;
                            closest_sphere_index = i;
                            closest_intersection_point = intersection_point;
                        }
                    }
                }

                //set the pixel color of the screen based on the color of the closest sphere if there is one
                if (closest_sphere_index != -1)
                {
                    //add the color of the closest sphere to the pixel color
                    add_vectors(&pixel_color, &(scene->spheres[closest_sphere_index].material.color));
                    bounces++;

                    //compute the new ray direction
                    Vector normal;
                    set_vector(&normal,
                               closest_intersection_point.x - scene->spheres[closest_sphere_index].center.x,
                               closest_intersection_point.y - scene->spheres[closest_sphere_index].center.y,
                               closest_intersection_point.z - scene->spheres[closest_sphere_index].center.z);
                    normalize_vector(&normal);
                    reflect_vector(&ray.direction, &normal);
                    normalize_vector(&ray.direction);
                    set_vector((Vector *)&ray.origin, closest_intersection_point.x, closest_intersection_point.y, closest_intersection_point.z);

                    //push the starting point of the ray outside the surface of the sphere
                    Vector offset = scale_vector_copy(&ray.direction, 0.001);
                    add_vectors((Vector *)&ray.origin, &offset);
                }
                else
                {
                    //add sky color to the pixel
                    add_vectors(&pixel_color, &SKY_COLOR);
                    still_going = 0;
                }
            }

            //set the color of this pixel in the screen
            // if (bounces > 0)
            {
                scale_vector(&pixel_color, 1.0 / (bounces + 1));
            }
            screen->pixels[row * screen->width + column] = pixel_color;
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

// //function to get the current linux epoch time in milliseconds
// long get_time_millis()
// {
//     struct timespec spec;
//     clock_gettime(CLOCK_REALTIME, &spec);
//     return spec.tv_sec * 1000 + spec.tv_nsec / 1000000;
// }

//test/create a simple scene with a single sphere in the middle and the camera looking directly at it
int main()
{
    //seed the random number generator
    srand(time(0));

    //get the time that the program starts
    struct timespec start;
    timespec_get(&start, TIME_UTC);

//create a list of 6 spheres 1 for each direction in 3D
#define NUM_SPHERES 6
    Sphere spheres[NUM_SPHERES] = {
        {.center = {.x = 0.25, .y = 0.0, .z = 0.0}, .material = {.color = {.x = 1.0, .y = 1.0, .z = 1.0}}, .radius = 0.125},
        {.center = {.x = 0.0, .y = 0.25, .z = 0.0}, .material = {.color = {.x = 0.5, .y = 0.5, .z = 0.5}}, .radius = 0.125},
        {.center = {.x = 0.0, .y = 0.0, .z = 0.25}, .material = {.color = {.x = 0.0, .y = 0.0, .z = 0.0}}, .radius = 0.125},
        {.center = {.x = -0.25, .y = 0.0, .z = 0.0}, .material = {.color = {.x = 0.0, .y = 1.0, .z = 1.0}}, .radius = 0.125},
        {.center = {.x = 0.0, .y = -0.25, .z = 0.0}, .material = {.color = {.x = 1.0, .y = 0.0, .z = 1.0}}, .radius = 0.125},
        {.center = {.x = 0.0, .y = 0.0, .z = -0.25}, .material = {.color = {.x = 1.0, .y = 1.0, .z = 0.0}}, .radius = 0.125},
    };

    //create a camera looking at the sphere from 2 meters away
    //camera position will be set inside the while loop
    Camera camera;
    init_camera(&camera);

    //create a screen of size SCREEN_WIDTH x SCREEN_HEIGHT with statically allocated pixels
    Vector pixels[SCREEN_HEIGHT * SCREEN_WIDTH];
    Screen screen = {.pixels = (Vector *)pixels, .width = SCREEN_WIDTH, .height = SCREEN_HEIGHT};

    //loop forever. set each pixel on the screen to a random color and draw the screen to the terminal
    struct timespec ts; //keep track of the system time for computing the duration of frames
    // double t = 0.0; //number of seconds elapsed since program started
    while (1)
    {
        //get the current system time, get the number of nanoseconds elapsed since the start of the program
        timespec_get(&ts, TIME_UTC);
        long long start_nanos = (ts.tv_sec - start.tv_sec) * 1000000000 + ts.tv_nsec - start.tv_nsec;

        //seconds that this frame occurred at
        double t = (double)start_nanos / 1000000000.0;

        //set the screen to the background color
        //set each pixel on the screen to a random color
        for (int i = 0; i < screen.height; i++)
        {
            for (int j = 0; j < screen.width; j++)
            {
                set_vector(&(screen.pixels[i * screen.width + j]), random_number(), random_number(), random_number());
            }
        }
        set_screen_color(&screen, BACKGROUND_COLOR);

        //project the scene onto the screen
        project_scene(&(Scene){.camera = camera, .spheres = spheres, .num_spheres = NUM_SPHERES}, &screen);

        //draw the screen to the terminal
        draw_screen(&screen);

        //construct the transform of the camera. the camera should orbit the center of the scene
        Frame tf0, tf1;
        init_frame(&tf0);
        init_frame(&tf1);
        init_frame(&(camera.frame));
        rotate_basis_y(&tf0.basis, 2.0 * PI * t * 0.3);
        rotate_basis_x(&tf0.basis, 2.0 * PI * t * 0.5);
        Vector root_to_camera = {.x = 0.0, .y = 0.0, .z = 5.0};
        add_vectors((Vector *)&tf1.origin, &root_to_camera);
        transform_frame(&camera.frame, &tf1);
        transform_frame(&camera.frame, &tf0);

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

        //print the number of nanoseconds at the bottom of the terminal, and then move the cursor back to the top
        // printf("%lld\n", program_nanos);
        // printf("\033[1;1H");
    }
}

// //struct for a light
// typedef struct
// {
//     Point position;
//     Vector color;
// } Light;

// //struct for a camera
// typedef struct
// {

// } Camera;

// //struct for a material
// typedef struct
// {
//     Vector color;
//     double ambient;
//     double diffuse;
//     double specular;
//     double shininess;
// } Material;

// //struct for a scene
// typedef struct
// {
//     Sphere *spheres;
//     int num_spheres;
//     Light *lights;
//     int num_lights;
//     Material *materials;
//     int num_materials;
// } Scene;

//struct for a camera
// typedef struct
// {
//     Point position;
//     Vector direction;
//     Vector up;
//     double fov;
// } Camera;

// //struct for a pixel
// typedef struct
// {
//     int x;
//     int y;
//     Vector color;
// } Pixel;

// //struct for a screen
// typedef struct
// {
//     Pixel *pixels;
//     int width;
//     int height;
// } Screen;

//struct for a ray intersection
// typedef struct
// {
//     double t;
//     Vector normal;
//     Material material;
// } Intersection;

//function prototypes
// void init_screen(Screen *screen);
// void init_scene(Scene *scene);
// void init_camera(Camera *camera);
// void init_light(Light *light);
// void init_material(Material *material);
// void init_sphere(Sphere *sphere);
// void init_ray(Ray *ray, Point origin, Vector direction);
// void init_intersection(Intersection *intersection);
// void init_point(Point *point, double x, double y, double z);
// void init_vector(Vector *vector, double x, double y, double z);
// void init_pixel(Pixel *pixel, int x, int y, Vector color);

// //function to print a pixel to the screen
// void print_pixel(Pixel *pixel, int color_code)
// {
//     printf("\033[%d;%dH", pixel->y, pixel->x);
//     printf("\033[48;5;%dm", color_code);
//     printf(" ");
//     printf("\033[0m");
// }

// //function to print a screen to the terminal
// void print_screen(Screen *screen)
// {
//     int i, j;
//     for (i = 0; i < screen->height; i++)
//     {
//         for (j = 0; j < screen->width; j++)
//         {
//             print_pixel(&screen->pixels[i * screen->width + j], screen->pixels[i * screen->width + j].color.x);
//         }
//         printf("\n");
//     }
// }

// //function to create a scene
// void create_scene(Scene *scene)
// {
//     //create the scene
//     init_scene(scene);

//     //create the camera
//     Camera camera;
//     init_camera(&camera);

//     //create the light
//     Light light;
//     init_light(&light);

//     //create the material
//     Material material;
//     init_material(&material);

//     //create the sphere
//     Sphere sphere;
//     init_sphere(&sphere);

//     //create the ray
//     Ray ray;
//     init_ray(&ray, camera.position, camera.direction);

//     //create the intersection
//     Intersection intersection;
//     init_intersection(&intersection);

//     //create the point
//     Point point;
//     init_point(&point, 0, 0, 0);

//     //create the vector
//     Vector vector;
//     init_vector(&vector, 0, 0, 0);

//     //create the pixel
//     Pixel pixel;
//     init_pixel(&pixel, 0, 0, vector);

//     //create the screen
//     Screen screen;
//     init_screen(&screen);

//     //create the scene
//     scene->num_spheres = 1;
//     scene->spheres = (Sphere *)malloc(sizeof(Sphere) * scene->num_spheres);
//     scene->spheres[0] = sphere;

//     scene->num_lights = 1;
//     scene->lights = (Light *)malloc(sizeof(Light) * scene->num_lights);
//     scene->lights[0] = light;

//     scene->num_materials = 1;
//     scene->materials = (Material *)malloc(sizeof(Material) * scene->num_materials);
//     scene->materials[0] = material;
// }

// //function to create a camera
// void init_camera(Camera *camera)
// {
//     //create the camera
//     camera->position.x = 0;
//     camera->position.y = 0;
//     camera->position.z = 0;

//     camera->direction.x = 0;
//     camera->direction.y = 0;
//     camera->direction.z = 1;

//     camera->up.x = 0;
//     camera->up.y = 1;
//     camera->up.z = 0;

//     camera->fov = PI / 3;
// }

// //function to create a light
// void init_light(Light *light)
// {
//     //create the light
//     light->position.x = 0;
//     light->position.y = 0;
//     light->position.z = 0;

//     light->color.x = 1;
//     light->color.y = 1;
//     light->color.z = 1;
// }

// //function to create a material
// void init_material(Material *material)
// {
//     //create the material
//     material->color.x = 1;
//     material->color.y = 1;
//     material->color.z = 1;

//     material->ambient = 0.1;
//     material->diffuse = 0.9;
//     material->specular = 0.9;
//     material->shininess = 200;
// }

// //function to create a sphere
// void init_sphere(Sphere *sphere)
// {
//     //create the sphere
//     sphere->center.x = 0;
//     sphere->center.y = 0;
//     sphere->center.z = 0;

//     sphere->radius = 1;

//     sphere->color.x = 1;
//     sphere->color.y = 1;
//     sphere->color.z = 1;
// }

// //function to create a ray
// void init_ray(Ray *ray, Point origin, Vector direction)
// {
//     //create the ray
//     ray->origin = origin;
//     ray->direction = direction;
// }

// //function to create an intersection
// void init_intersection(Intersection *intersection)
// {
//     //create the intersection
//     intersection->t = 0;
//     intersection->normal.x = 0;
//     intersection->normal.y = 0;
//     intersection->normal.z = 0;
//     intersection->material.color.x = 0;
//     intersection->material.color.y = 0;
//     intersection->material.color.z = 0;
//     intersection->material.ambient = 0;
//     intersection->material.diffuse = 0;
//     intersection->material.specular = 0;
//     intersection->material.shininess = 0;
// }

// //function to create a point
// void init_point(Point *point, double x, double y, double z)
// {
//     //create the point
//     point->x = x;
//     point->y = y;
//     point->z = z;
// }

// //function to create a vector
// void init_vector(Vector *vector, double x, double y, double z)
// {
//     //create the vector
//     vector->x = x;
//     vector->y = y;
//     vector->z = z;
// }

// //function to create a pixel
// void init_pixel(Pixel *pixel, int x, int y, Vector color)
// {
//     //create the pixel
//     pixel->x = x;
//     pixel->y = y;
//     pixel->color = color;
// }

// //function to create a screen
// void init_screen(Screen *screen)
// {
//     //create the screen
//     screen->width = 80;
//     screen->height = 25;
//     screen->pixels = (Pixel *)malloc(sizeof(Pixel) * screen->width * screen->height);
// }

// //function to create a scene
// void init_scene(Scene *scene)
// {
//     //create the scene
//     scene->num_spheres = 0;
//     scene->spheres = NULL;

//     scene->num_lights = 0;
//     scene->lights = NULL;

//     scene->num_materials = 0;
//     scene->materials = NULL;
// }

// //function to create a ray
// void create_ray(Ray *ray, Point origin, Vector direction)
// {
//     //create the ray
//     ray->origin = origin;
//     ray->direction = direction;
// }

// //function to create a sphere
// void create_sphere(Sphere *sphere, Point center, double radius, Vector color)
// {
//     //create the sphere
//     sphere->center = center;
//     sphere->radius = radius;
//     sphere->color = color;
// }

// //function to create a material
// void create_material(Material *material, Vector color, double ambient, double diffuse, double specular, double shininess)
// {
//     //create the material
//     material->color = color;
//     material->ambient = ambient;
//     material->diffuse = diffuse;
//     material->specular = specular;
//     material->shininess = shininess;
// }

// //function to create a light
// void create_light(Light *light, Point position, Vector color)
// {
//     //create the light
//     light->position = position;
//     light->color = color;
// }

// //function to create an intersection
// void create_intersection(Intersection *intersection, double t, Vector normal, Material material)
// {
//     //create the intersection
//     intersection->t = t;
//     intersection->normal = normal;
//     intersection->material = material;
// }

// //function to create a point
// void create_point(Point *point, double x, double y, double z)
// {
//     //create the point
//     point->x = x;
//     point->y = y;
//     point->z = z;
// }

// //function to create a vector
// void create_vector(Vector *vector, double x, double y, double z)
// {
//     //create the vector
//     vector->x = x;
//     vector->y = y;
//     vector->z = z;
// }

// //function to create a pixel
// void create_pixel(Pixel *pixel, int x, int y, Vector color)
// {
//     //create the pixel
//     pixel->x = x;
//     pixel->y = y;
//     pixel->color = color;
// }

// //function to create a screen
// void create_screen(Screen *screen)
// {
//     //create the screen
//     screen->width = 80;
//     screen->height = 25;
//     screen->pixels = (Pixel *)malloc(sizeof(Pixel) * screen->width * screen->height);
// }

// //function to run the ray tracer
// void run_ray_tracer(Scene *scene, Camera *camera, Light *light, Screen *screen)
// {
//     //draw the scene to the terminal
//     draw_scene(scene, camera, light, screen);

//     //update the camera position by rotating the camera around the scene
//     update_camera(scene, camera);
// }