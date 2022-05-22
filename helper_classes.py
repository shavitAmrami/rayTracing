import numpy as np


# This function gets a vector and returns its normalized form.
def normalize(vector):
    return vector / np.linalg.norm(vector)


def check_if_shadowed(light_source, light_ray, objects):
    shadow_coefficient = 1
    min_distance, _ = light_ray.nearest_intersected_object(objects=objects)
    if isinstance(light_source, DirectionalLight):
        if 0 < min_distance < np.inf:
            shadow_coefficient = 0
    else:
    # if isinstance(light_source, PointlLight) or isinstance(light_source, SpotlLigh
        if 0 < min_distance < 1:
            shadow_coefficient = 0
    return shadow_coefficient


def calc_diffuse_component(light_ray, light_intensity, object_normal, nearest_object):
    diffuse_coefficient = nearest_object.diffuse
    surface_normal = object_normal
    dot_product = np.dot(surface_normal, normalize(light_ray.direction))
    diffuse_component = np.multiply(diffuse_coefficient, light_intensity) * dot_product
    return diffuse_component


def calc_specular_component(light_ray, light_intensity, object_normal, nearest_object, view_point_ray):
    specular_coefficient = nearest_object.specular
    shininess = nearest_object.shininess
    surface_normal = object_normal
    L_hat = reflected(vector=light_ray.direction, normal=surface_normal)
    scalar_calculation = np.dot(view_point_ray, L_hat) ** shininess
    specular_component = np.multiply(specular_coefficient, light_intensity) * scalar_calculation
    return specular_component


# This function gets a vector and the normal of the surface it hit
# This function returns the vector that reflects from the surface
def reflected(vector, normal):
    v = np.array([0, 0, 0])
    projection = np.dot(vector, normal) * normal
    v = np.subtract(vector, 2 * projection)
    return v


## Lights

class LightSource:

    def __init__(self, intensity):
        self.intensity = intensity


class DirectionalLight(LightSource):

    def __init__(self, intensity, direction):
        super().__init__(intensity)
        self.direction = direction

    # This function returns the ray that goes from the light source to a point
    def get_light_ray(self, intersection):
        return Ray(origin=intersection, direction=normalize(self.direction))

    # This function returns the distance from a point to the light source
    def get_distance_from_light(self, intersection):
        return np.inf

    # This function returns the light intensity at a point
    def get_intensity(self, intersection):
        return self.intensity


class PointLight(LightSource):

    def __init__(self, intensity, position, kc, kl, kq):
        super().__init__(intensity)
        self.position = np.array(position)
        self.kc = kc
        self.kl = kl
        self.kq = kq

    # This function returns the ray that goes from the light source to a point
    def get_light_ray(self, intersection):
        return Ray(origin=intersection, direction=normalize(self.position - intersection))

    # This function returns the distance from a point to the light source
    def get_distance_from_light(self, intersection):
        return np.linalg.norm(intersection - self.position)

    # This function returns the light intensity at a point
    def get_intensity(self, intersection):
        d = self.get_distance_from_light(intersection)
        return self.intensity / (self.kc + self.kl*d + self.kq * (d**2))


class SpotLight(LightSource):

    def __init__(self, intensity, position, direction, kc, kl, kq):
        super().__init__(intensity)
        self.position = np.array(position)
        self.central_direction = (-1) * normalize(np.array(direction))
        self.kc = kc
        self.kl = kl
        self.kq = kq

    # This function returns the ray that goes from the light source to a point
    def get_light_ray(self, intersection):
        return Ray(origin=intersection, direction=normalize(self.position - intersection))

    def get_distance_from_light(self, intersection):
        return np.linalg.norm(intersection - self.position)

    def get_intensity(self, intersection):
        d = self.get_distance_from_light(intersection)
        v = normalize(intersection - self.position)
        v_dot_central_direction = np.dot(v, self.central_direction)
        return (self.intensity * v_dot_central_direction) / (self.kc + self.kl*d + self.kq * (d**2))


class Ray:
    def __init__(self, origin, direction):
        self.origin = origin
        self.direction = direction

    # The function is getting the collection of objects in the scene and looks for the one with minimum distance.
    # The function should return the nearest object and its distance (in two different arguments)
    def nearest_intersected_object(self, objects):
        nearest_object = None
        min_distance = np.inf

        for scene_object in objects:
            t, optional_nearest_object = scene_object.intersect(self)
            if t is not None and t < min_distance and optional_nearest_object is not None:
            # if t < min_distance and optional_nearest_object is not None:
                nearest_object = optional_nearest_object
                min_distance = t
        return min_distance, nearest_object


class Object3D:

    def set_material(self, ambient, diffuse, specular, shininess, reflection):
        self.ambient = ambient
        self.diffuse = diffuse
        self.specular = specular
        self.shininess = shininess
        self.reflection = reflection


class Plane(Object3D):
    def __init__(self, normal, point):
        self.normal = np.array(normal)
        self.point = np.array(point)

    def intersect(self, ray: Ray):
        v = self.point - ray.origin
        t = (np.dot(v, self.normal) / np.dot(self.normal, ray.direction))
        if t > 0:
            return t, self
        else:
            return -1, None


class Triangle(Object3D):
    # Triangle gets 3 points as arguments
    def __init__(self, a, b, c):
        self.a = np.array(a)
        self.b = np.array(b)
        self.c = np.array(c)
        self.normal = self.compute_normal()

    def compute_normal(self):
        n = np.cross((self.b - self.a), (self.c - self.a))
        return n

    # Hint: First find the intersection on the plane
    # Later, find if the point is in the triangle using barycentric coordinates
    def intersect(self, ray: Ray):
        # Find the intersection point (P) with the plane that contains the triangle:
        triangle_plane_normal = self.normal
        triangle_plane = Plane(normal=triangle_plane_normal, point=self.a)
        t, _ = triangle_plane.intersect(ray=ray)

        if t == -1:  # Means that there is no intersection with the triangle plane
            return -1, None
        else:
            q = ray.origin + (t * ray.direction)
            first_condition = np.dot(np.cross((self.b - self.a), (q - self.a)), self.normal)
            second_condition = np.dot(np.cross((self.c - self.b), (q - self.b)), self.normal)
            third_condition = np.dot(np.cross((self.a - self.c), (q - self.c)), self.normal)

            if first_condition >= 0 and second_condition >= 0 and third_condition >= 0:
                return t, self
        return -1, None


class Sphere(Object3D):
    def __init__(self, center, radius: float):
        self.center = center
        self.radius = radius

    def intersect(self, ray: Ray):
        b = 2 * np.dot(ray.direction, ray.origin - self.center)
        c = np.linalg.norm(ray.origin - self.center) ** 2 - self.radius ** 2
        delta = b ** 2 - 4 * c
        if delta > 0:
            t1 = (-b + np.sqrt(delta)) / 2
            t2 = (-b - np.sqrt(delta)) / 2
            if t1 > 0 and t2 > 0:
                return min(t1, t2), self
        return -1, None

    def compute_normal_at_intersection(self, intersection):
        n = intersection - self.center
        return n


class Mesh(Object3D):
    # Mesh are defined by a list of vertices, and a list of faces.
    # The faces are triplets of vertices by their index number.
    def __init__(self, v_list, f_list):
        self.v_list = v_list
        self.f_list = f_list
        self.triangle_list = self.create_triangle_list()

    def create_triangle_list(self):
        l = []
        for face in self.f_list:
            a = self.v_list[face[0]]
            b = self.v_list[face[1]]
            c = self.v_list[face[2]]
            l.append(Triangle(a=a, b=b, c=c))
        return l

    def apply_materials_to_triangles(self):
        for t in self.triangle_list:
            t.set_material(self.ambient, self.diffuse, self.specular, self.shininess, self.reflection)

    # Hint: Intersect returns both distance and nearest object.
    # Keep track of both.
    def intersect(self, ray: Ray):
        intersection_list = []
        for triangle in self.triangle_list:
            curr_t, curr_triangle = triangle.intersect(ray=ray)
            intersection_list.append((curr_t, curr_triangle))

        # return min(intersection_list, key=lambda t: t[0])
        return ray.nearest_intersected_object(objects=self.triangle_list)

