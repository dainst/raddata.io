import geopandas as gpd
import numpy


class Shp2Rast:
    def __init__(self, fname: str):
        self.all_data = gpd.read_file(fname)
        try:
            self.geo = self.all_data['geometry']
        except KeyError:
            print('no gepmetry column present!')
        self.x, self.y = Shp2Rast.grab_coordsxy(self.geo)
        self.x_res = numpy.abs(numpy.nanmedian(self.x - numpy.roll(x, 1)))
        self.y_res = numpy.abs(numpy.nanmedian(self.y - numpy.roll(y, 1)))
        self.xmax = self.x.max()
        self.ymax = self.y.max()
        self.ymin = self.y.min()
        self.xmin = self.x.min()
        self.xi = numpy.arange(self.ymin, self.ymax, self.y_res)
        self.yi = numpy.arange(self.xmin, self.xmax, self.x_res)

    def interpolate(self, columnname: str):
        self.z = self.all_data.get(columnname).to_numpy()
        self.zi = Shp2Rast.scp_gdd_interp(self.x, self.y, self.z, self.xi, self.yi)

    @staticmethod
    def grab_coordsxy(coordinates: gpd.geoseries):
        x = []
        y = []
        for j in coordinates:
            if j.centroid[0] == j.centroid[2] and j.centroid[1] == j.centroid[3]:
                x.append(j.centroid[0])
                y.append(j.centroid[1])
            else:
                print('Uups')
        return numpy.asarray(x), numpy.asarray(y)

    @staticmethod
    def scp_gdd_interp(x: numpy.ndarray, y: numpy.ndarray, z: numpy.ndarray, xi: numpy.ndarray, yi: numpy.ndarray,
                       method="cubic") -> numpy.ndarray:
        """

        :param x:1d Vector of x values of the data points
        :param y:1d Vector of y values of the data points
        :param z: Values of the points
        :param xi:1d Vector of X values for interpolation
        :param yi: 1d Vector of Y values for interpolation
        :param method:
        :return: zint
        """
        grid_x, grid_y = numpy.meshgrid(xi, yi)  # interpolate at these points
        to_interp = numpy.c_[grid_x.ravel(), grid_y.ravel()]
        points = numpy.c_[x.ravel(), y.ravel()]
        if method == "cubic":
            zint = griddata(points, z.ravel(), to_interp, method='cubic')
        if method == "nearest":
            zint = griddata(points, z.ravel(), to_interp, method='nearest')
        if method == "linear":
            zint = griddata(points, z.ravel, to_interp, method='linear')
        zint = zint.reshape(grid_x.shape)
        return grid_x, grid_y, zint
