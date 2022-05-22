from helper_classes import *
import matplotlib.pyplot as plt


def render_scene(camera, ambient, lights, objects, screen_size, max_depth):
    width, height = screen_size
    ratio = float(width) / height
    screen = (-1, 1 / ratio, 1, -1 / ratio)  # left, top, right, bottom

    image = np.zeros((height, width, 3))

    for i, y in enumerate(np.linspace(screen[1], screen[3], height)):
        for j, x in enumerate(np.linspace(screen[0], screen[2], width)):
            pixel = np.array([x, y, 0])
            color = np.zeros(3)
            origin = camera
            ray_direction = normalize(pixel - origin)
            reflection = 1

            for level in range(max_depth):
                reflectance = np.zeros(3)
                ray = Ray(origin, ray_direction)
                min_t, nearest_object = ray.nearest_intersected_object(objects=objects)
                if nearest_object is None:  # There is no intersection -> pixel's color will be black
                    break

                intersection = origin + min_t * ray_direction
                if isinstance(nearest_object, Sphere):
                    object_normal = normalize(nearest_object.compute_normal_at_intersection(intersection))
                else:
                    object_normal = normalize(nearest_object.normal)
                shifted_intersection_point = intersection + 1e-6 * object_normal

                ambient_component = np.multiply(ambient, nearest_object.ambient)
                total_diffuse_component = np.zeros(3)
                total_specular_component = np.zeros(3)
                view_point_ray = normalize(origin - intersection)

                for light_source in lights:
                    light_ray = light_source.get_light_ray(shifted_intersection_point)
                    light_intensity = light_source.get_intensity(shifted_intersection_point)
                    shadow_coefficient = check_if_shadowed(light_source, light_ray, objects)

                    # calculate diffuse component
                    diffuse_by_current_light_source = calc_diffuse_component(light_ray, light_intensity,
                                                                             object_normal, nearest_object)
                    total_diffuse_component = (np.add(total_diffuse_component,
                                                      diffuse_by_current_light_source)) * shadow_coefficient

                    # calc specular component
                    specular_by_current_light_source = calc_specular_component(light_ray, light_intensity,
                                                                               object_normal, nearest_object,
                                                                               view_point_ray)
                    total_specular_component = (np.add(total_specular_component,
                                                       specular_by_current_light_source)) * shadow_coefficient

                    reflectance = np.add(reflectance, total_diffuse_component + total_specular_component)
                reflectance += ambient_component
                color = np.add(color, reflection * reflectance)
                reflection *= nearest_object.reflection

                origin = shifted_intersection_point
                ray_direction = reflected(vector=ray_direction, normal=object_normal)

            # We clip the values between 0 and 1 so all pixel values will make sense.
            image[i, j] = np.clip(color, 0, 1)

    return image


# Write your own objects and lights
def your_own_scene():
    camera = np.array([0, 0, 1])

    # Lights:
    light_a = DirectionalLight(intensity=np.array([1, 0.2, 1]), direction=np.array([1, -1, 0.9]))
    light_b = SpotLight(intensity=np.array([1, 1, 0]), position=np.array([0, 0.5, 0]), direction=([0, 0.5, 1]),
                        kc=0.3, kl=0.7, kq=0.1)
    lights = [light_a, light_b]
    # ambient = np.array([0.1, 0.1, 0.1])

    # Objects:
    plane = Plane([0, 0, 1], [0, 0, -3])
    plane.set_material([0, 0.5, 2], [0, 1, 1], [1, 1, 1], 200, 0)
    triangle = Triangle([1, 0.5, -2], [0, 1, -1.5], [0, -1, -1])
    triangle.set_material([1, 0, 0], [1, 0.753, 0.78], [0, 0, 0], 100, 0.5)
    sphere = Sphere([-0.5, -0.3, -0.1], 0.2)
    sphere.set_material([0, 1, 0], [1, 1, 0], [0.3, 0.3, 0.3], 100, 0)
    objects = [plane, triangle, sphere]

    return camera, lights, objects
