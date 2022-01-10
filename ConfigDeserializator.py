import xml.etree.ElementTree as ET


def read_config(config):
    tree = ET.parse(config)
    return tree.getroot()


def read_node(tree: ET.Element, node):
    return tree.find(node)


def read_dem_path(area: ET.Element):
    dem = area.find("DemImage")
    return dem.find("PathToImage").text


def read_radar_image_attrs(radar_node: ET.Element):
    area_attrs = {'subswatch': radar_node.find("Subswatch").text,
                  'first_burst_index': radar_node.find("FirstBurstIndex").text,
                  'last_burst_index': radar_node.find("LastBurstIndex").text
                  }
    return area_attrs


def read_images_path(node: ET.Element):
    images_path = {'flood_image': node.find("PathToFloodImage").text,
                   'non_flood_image': node.find("PathToNonFloodImage").text
                   }
    return images_path
