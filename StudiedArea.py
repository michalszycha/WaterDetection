import ConfigDeserializator


class Images:
    def __init__(self, tree, node_str):
        self._node = ConfigDeserializator.read_node(tree, node_str)
        self._flood_image_path = ConfigDeserializator.read_images_path(self._node)['flood_image']
        self._non_flood_image_path = ConfigDeserializator.read_images_path(self._node)['non_flood_image']

    @property
    def flood_image_path(self):
        return self._flood_image_path

    @property
    def non_flood_image_patch(self):
        return self._non_flood_image_path


class StudiedArea:
    def __init__(self, config):
        self._config = config
        self._root = ConfigDeserializator.read_config(self._config)
        self._area = ConfigDeserializator.read_node(self._root, "area")
        self._area_name = self._area.get("name")
        self._area_coordinates = ConfigDeserializator.read_node(self._area, "AreaOfInterest").text
        self._dem_image = self.DemImage(self._area)
        self._radar_images = self.RadarImages(self._area)
        self._optical_images = self.OpticalImages(self._area, "OpticalImages")

    @property
    def area_name(self):
        return self._area_name

    @property
    def area_coordinates(self):
        return self._area_coordinates

    @property
    def dem_image(self):
        return self._dem_image

    @property
    def radar_images(self):
        return self._radar_images

    @property
    def optical_images(self):
        return self._optical_images

    class DemImage:
        def __init__(self, area):
            self._path_to_image = ConfigDeserializator.read_dem_path(area)

        @property
        def path_to_image(self):
            return self._path_to_image

    class RadarImages:
        def __init__(self, area):
            self._radar_node = ConfigDeserializator.read_node(area, "RadarImages")
            self._image_attrs = ConfigDeserializator.read_radar_image_attrs(self._radar_node)
            self._subswatch = self._image_attrs['subswatch']
            self._first_burst_index = self._image_attrs['first_burst_index']
            self._last_burst_index = self._image_attrs['last_burst_index']
            self._slc = self.Slc(self._radar_node, "Slc")
            self._grd = self.Grd(self._radar_node, "Grd")

        @property
        def subswatch(self):
            return self._subswatch

        @property
        def first_burst_index(self):
            return self._first_burst_index

        @property
        def last_burst_index(self):
            return self._last_burst_index

        @property
        def slc(self):
            return self._slc

        @property
        def grd(self):
            return self._grd

        class Slc(Images):
            pass

        class Grd(Images):
            pass

    class OpticalImages(Images):
        pass
