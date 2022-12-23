#include <SFML/Graphics.hpp>
#include <queue>
#include <fstream>

using namespace std;
using namespace sf;

class MyPoint;
class MySegment;
class MyRectangle;
class Scene;
class Graph;

class MyPoint
{
public:
	MyPoint(double x, double y) : _x(x), _y(y) {}

	double get_x() const { return _x; }
	double get_y() const { return _y; }

	void set_x(double x) { _x = x; }
	void set_y(double y) { _y = y; }

	double distance(const MyPoint& other) const
	{
		return sqrt(pow(other.get_x() - _x, 2) + pow(other.get_y() - _y, 2));
	}

	unique_ptr<MyPoint> find_closest(vector<MyPoint> points)
	{
		unique_ptr<MyPoint> result = NULL;
		double min = numeric_limits<int>::max();
		for (auto& p : points)
			if (this->distance(p) < min) {
				min = this->distance(p);
				result = make_unique<MyPoint>(p);
			}

		return result;
	}

	unique_ptr<MyRectangle> find_closest_rectangle(vector<MyRectangle> rectangles);

	bool in(MyRectangle& rectangle);

	bool operator<(const MyPoint& p) const
	{
		return (_x < p.get_x()) || (equal(_x, p.get_x()) && _y < p.get_y());
	};

	bool operator==(const MyPoint& p) const
	{
		return equal(_x, p.get_x()) && equal(_y, p.get_y());
	};

	friend ofstream& operator<<(ofstream& file, const MyPoint& p)
	{
		file << p._x << ' ' << p._y;

		return file;
	};

private:
	double _x;
	double _y;

	bool equal(double a, double b) const
	{
		static const double EPS = 1e-8;
		return abs(a - b) < EPS;
	}
};

class MySegment
{
public:
	MySegment(MyPoint p1, MyPoint p2) : _p1(p1), _p2(p2) {}

	MyPoint getA() const { return _p1; }

	MyPoint getB() const { return _p2; }

	double dot(const MySegment& s) const
	{
		pair<double, double> u = this->to_vector();
		pair<double, double> v = s.to_vector();
		return u.first * v.first + u.second * v.second;
	}

	double cross(const MySegment& s) const
	{
		pair<double, double> u = this->to_vector();
		pair<double, double> v = s.to_vector();
		return u.first * v.second - u.second * v.first;
	}

	pair<double, double> to_vector() const
	{
		return make_pair(_p2.get_x() - _p1.get_x(), _p2.get_y() - _p1.get_y());
	}

private:
	MyPoint _p1;
	MyPoint _p2;
};

class MyRectangle
{
public:
	MyRectangle(MyPoint p1, MyPoint p2) : _p1(p1), _p2(p2)
	{
		if (_p1.get_x() > _p2.get_x())
			swap(_p1, _p2);
		if (_p1.get_y() > _p2.get_y())
		{
			double ty = _p1.get_y();
			_p1 = MyPoint(_p1.get_x(), _p2.get_y());
			_p2 = MyPoint(_p2.get_x(), ty);
		}
	}

	double getXsize() const
	{
		return _p2.get_x() - _p1.get_x() + 1;
	}

	double getYsize() const
	{
		return _p2.get_y() - _p1.get_y() + 1;
	}

	MyPoint getMainPoint() const
	{
		return _p1;
	}

	MyPoint getCentralPoint() const
	{
		return MyPoint((_p1.get_x() + _p2.get_x()) / 2, (_p1.get_y() + _p2.get_y()) / 2);
	}

	vector<MyPoint> getVertices(float margin = 0) const
	{
		vector<MyPoint> vertices;

		vertices.push_back(MyPoint(_p1.get_x() - margin - 1, _p1.get_y() - margin - 1));
		vertices.push_back(MyPoint(_p2.get_x() + margin + 1, _p1.get_y() - margin - 1));
		vertices.push_back(MyPoint(_p2.get_x() + margin + 1, _p2.get_y() + margin + 1));
		vertices.push_back(MyPoint(_p1.get_x() - margin - 1, _p2.get_y() + margin + 1));

		return vertices;
	};

	bool operator==(const MyRectangle& r) const
	{
		return _p1 == r.getMainPoint() && this->getXsize() == r.getXsize() && this->getYsize() == r.getYsize();
	};

	friend ofstream& operator<<(ofstream& file, const MyRectangle& r)
	{
		file << r._p1 << ' ';
		file << r._p2;

		return file;
	};

private:
	MyPoint _p1;
	MyPoint _p2;
};

unique_ptr<MyRectangle> MyPoint::find_closest_rectangle(vector<MyRectangle> rectangles)
{
	unique_ptr<MyRectangle> result = NULL;
	double min = numeric_limits<int>::max();
	for (auto& p : rectangles)
	{
		double distance = this->distance(p.getCentralPoint());
		if (distance < min) {
			min = distance;
			result = make_unique<MyRectangle>(p);
		}
	}

	return result;
}

bool MyPoint::in(MyRectangle& rectangle)
{
	int pos = 0;
	if (_x < rectangle.getMainPoint().get_x())	pos += 1;
	if (_x > rectangle.getMainPoint().get_x() + rectangle.getXsize())	pos += 1 << 1;
	if (_y < rectangle.getMainPoint().get_y())	pos += 1 << 2;
	if (_y > rectangle.getMainPoint().get_y() + rectangle.getYsize())	pos += 1 << 3;
	if (!pos)
		return true;
	return false;
}

class Scene
{
public:
	Scene()
	{
		if (!font.loadFromFile("resources/arial.ttf"))
			throw invalid_argument("Can't find text font 'arial.ttf'\n");

		edit_mode_indicator.setFont(font);
		edit_mode_indicator.setCharacterSize(16);
		edit_mode_indicator.setPosition(Vector2f(6.0f, 27.0f));
		edit_mode_indicator.setOutlineColor(Color(255, 51, 0));
		edit_mode_indicator.setOutlineThickness(2.f);
		edit_mode_indicator.setFillColor(Color(255, 255, 0));
		edit_mode_indicator.setString("BUILDING  MODE  ACTIVE");

		preview_rect.setFillColor(Color(169, 169, 169, 100));

		help_underlay.setOutlineColor(Color::Black);
		help_underlay.setOutlineThickness(1.0f);
		help_underlay.setFillColor(Color::White);
		help_underlay.setPosition(Vector2f(1.f, 315.f));
		help_underlay.setSize(Vector2f(514.f, 404.f));

		help.setFont(font);
		help.setOutlineColor(Color::Black);
		help.setOutlineThickness(1.0f);
		help.setFillColor(Color::White);
		help.setCharacterSize(16);
		help.setPosition(5.0f, 320.0f);
		help.setString(
			"F1 - toggle help\n\n"
			"Shift + Mouse1 - set position of a robot\n"
			"Shift + Mouse2 - set destination\n"
			"Shift + MouseWheel - change size of a robot\n"
			"Mouse1 - add intermidiate\n"
			"Mouse2 - delete intermidiate\n"
			"Shift + C - delete all intermidiates\n"
			"H - toggle graph edges\n"
			"B - toggle building mode\n"
			"In building mode :\n"
			"\tMouse1 - set first vertice of a rectangle\n"
			"\tClick Mouse1 elsewhere or drag it to set second vertice of a rectangle\n"
			"\tMouse2 - cancel (first vertice is set)\n"
			"\tMouse2 - delete rectangle\n"
			"\tShift + C - delete all rectangles\n"
			"\tEscape - exit building mode\n"
			"Shift + D - clear the scene (make sure you saved it)\n\n"
			"F5 - save the scene\n"
			"F9 - load the scene\n"
			"Escape - exit"
		);

		path_lenght.setFont(font);
		path_lenght.setOutlineColor(Color::Black);
		path_lenght.setOutlineThickness(1.0f);
		path_lenght.setFillColor(Color::White);
		path_lenght.setCharacterSize(16);
		path_lenght.setPosition(5.0f, 5.0f);
		path_lenght.setString("Place robot and destination or load a save");
	}

	void add_rectangle(MyRectangle& rect)
	{
		RectangleShape rectangle(Vector2f(float(rect.getXsize()), float(rect.getYsize())));
		rectangle.setPosition(Vector2f(float(rect.getMainPoint().get_x()), float(rect.getMainPoint().get_y())));
		rectangle.setFillColor(Color(169, 169, 169));
		rectangles.push_back(rectangle);
	}

	void remove_rectangle(MyRectangle& r)
	{
		auto is_same = [r](RectangleShape& rect)
		{
			return rect.getPosition().x == r.getMainPoint().get_x() && rect.getPosition().y == r.getMainPoint().get_y()
				&& rect.getSize().x == r.getXsize() && rect.getSize().y == r.getYsize();
		};
		rectangles.erase(remove_if(rectangles.begin(), rectangles.end(), is_same), rectangles.end());
	}

	void upd_preview_rect()
	{
		preview_rect.setSize(Vector2f(0, 0));
		preview_rect.setPosition(Vector2f(0, 0));
	}

	void upd_preview_rect(shared_ptr<MyPoint> first, int x, int y)
	{
		preview_rect.setSize(Vector2f(
			float(abs(first->get_x() - x) + 1),
			float(abs(first->get_y() - y) + 1)
		));
		preview_rect.setPosition(Vector2f(
			float(first->get_x() > x ? x : first->get_x()),
			float(first->get_y() > y ? y : first->get_y())
		));
	}

	void set_robot() { robot = CircleShape(); }

	void set_robot(shared_ptr<MyPoint> p, const float radius = 5)
	{
		robot.setRadius(0.5f + radius);
		robot.setPosition(Vector2f(float(p->get_x() - int(radius + 0.5f)), float(p->get_y() - int(0.5f + radius))));
		robot.setFillColor(Color(31, 122, 31));
	}

	void set_destination() { destination = CircleShape(); }

	void set_destination(shared_ptr<MyPoint> p)
	{
		destination.setRadius(5.5);
		destination.setPosition(Vector2f(float(p->get_x() - 5), float(p->get_y() - 5)));
		destination.setFillColor(Color(179, 0, 0));
	}

	void add_dot(MyPoint& p)
	{
		CircleShape dot(5.5);
		dot.setPosition(Vector2f(float(p.get_x() - 5), float(p.get_y() - 5)));
		dot.setFillColor(Color(0, 102, 255));
		dots.push_back(dot);
	}

	void remove_dot(MyPoint& p)
	{	
		auto is_same = [p](CircleShape& dot) { return dot.getPosition().x + 5 == p.get_x() && dot.getPosition().y + 5 == p.get_y(); };
		dots.erase(remove_if(dots.begin(), dots.end(), is_same), dots.end());
	}

	void render_path(pair<pair<vector<MySegment>, vector<MySegment>>, float> graph_path_pathLenght)
	{
		graph_lines.clear();
		for (auto& edge : graph_path_pathLenght.first.first)
		{
			VertexArray new_line(Lines, 2);
			new_line[0].position = Vector2f(float(edge.getA().get_x() + 1), float(edge.getA().get_y()));
			new_line[1].position = Vector2f(float(edge.getB().get_x() + 1), float(edge.getB().get_y() + 1));
			new_line[0].color = Color(0, 179, 143, 70);
			new_line[1].color = Color(0, 179, 143, 70);
			graph_lines.push_back(new_line);
		}

		path_lines.clear();
		for (auto& edge : graph_path_pathLenght.first.second)
		{
			VertexArray new_line(Lines, 2);
			new_line[0].position = Vector2f(float(edge.getA().get_x() + 1), float(edge.getA().get_y()));
			new_line[1].position = Vector2f(float(edge.getB().get_x() + 1), float(edge.getB().get_y() + 1));
			new_line[0].color = Color(0, 102, 255);
			new_line[1].color = Color(0, 102, 255);
			path_lines.push_back(new_line);
		}

		if (graph_path_pathLenght.first.second.size())
			path_lenght.setString("Path lenght: " + to_string(graph_path_pathLenght.second));
		else
			path_lenght.setString("Path is impossible");

	}

	void clear()
	{
		rectangles.clear();
		robot = CircleShape();
		destination = CircleShape();
		dots.clear();
		graph_lines.clear();
		path_lines.clear();
		this->upd_preview_rect();
		path_lenght.setString("Place robot and destination or load a save");
	}

	vector<RectangleShape> rectangles;
	CircleShape robot;
	CircleShape destination;
	vector<CircleShape> dots;

	vector<VertexArray> graph_lines;
	vector<VertexArray> path_lines;

	RectangleShape preview_rect;
	Text edit_mode_indicator;

	Font font;
	Text help;
	RectangleShape help_underlay;
	Text path_lenght;
};

class Graph {
public:
	Graph() {}

	void add_vertex(int vertex) {
		if (!has_vertex(vertex)) {
			vertices[vertex] = std::map<int, double>();
		}
	}

	void add_edge(int start_vertex, int end_vertex, double weight) {
		add_vertex(start_vertex);
		add_vertex(end_vertex);
		vertices[start_vertex][end_vertex] = weight;
		vertices[end_vertex][start_vertex] = weight;
	}

	std::vector<int> get_vertices() const {
		std::vector<int> result;
		for (const auto& p : vertices) {
			result.push_back(p.first);
		}
		return result;
	}

	std::vector <pair<pair<int, int>, double>> get_edges() const {
		vector <pair<pair<int, int>, double>> result;
		for (const auto& u : vertices) {
			for (const auto& v : u.second) {
				auto is_vu = [u, v](pair<pair<int, int>, double> e) { return e.first.first == v.first && e.first.second == u.first; };
				if (std::find_if(result.begin(), result.end(), is_vu) == result.end())
					result.push_back(std::make_pair(std::make_pair(u.first, v.first), v.second));
			}
		}
		return result;
	}

	std::vector<int> get_adjacent_vertices(int src_vertex) const {
		const auto it = vertices.find(src_vertex);
		if (it == vertices.end()) {
			return std::vector<int> {};
		}
		vector<int> result;
		for (const auto& p : it->second) {
			result.push_back(p.first);
		}
		return result;
	}

	vector<pair<int, double>> get_adjacent_edges(int src_vertex) const {
		const auto it = vertices.find(src_vertex);
		if (it == vertices.end()) {
			return vector<pair<int, double>> {};
		}
		vector<pair<int, double>> result;
		for (const auto& p : it->second) {
			result.push_back(make_pair(p.first, p.second));
		}
		return result;
	}

	bool has_vertex(int vertex) const {
		return (vertices.find(vertex) != vertices.end());
	}

	bool has_edge(int start_vertex, int end_vertex) const {
		const auto it = vertices.find(start_vertex);
		if (it == vertices.end()) {
			return false;
		}
		return (it->second.find(end_vertex) != it->second.end());
	}

private:
	map<int, map<int, double>> vertices;
};

vector<int> build_shortest_path(vector<int> parent, int start_vertex, int end_vertex)
{
	vector<int> result;

	if (start_vertex != end_vertex)
	{
		result.push_back(end_vertex);

		for (int cur_v = end_vertex; cur_v != start_vertex;)
		{
			result.push_back(parent[cur_v]);
			cur_v = parent[cur_v];
		}
		reverse(result.begin(), result.end());
	}

	return result;
}

pair<vector<int>, double> shortest_path(const Graph& graph, int start_vertex, int end_vertex)
{
	// Return shortest path in the graph from start vertex to end vertex as array of vertices.
	// First item in the result should be start vertex, last - end vertex.
	// Return empty array if there is no path.

	vector<int> vertices = graph.get_vertices();

	if (vertices.size())
	{
		vector<double> distance(vertices.size(), double(numeric_limits<int>::max()));
		vector<int> parent(vertices.size());
		distance[start_vertex] = 0;

		priority_queue<pair<double, int>> queue;
		queue.push(make_pair(0, start_vertex));

		while (!queue.empty())
		{
			int u = queue.top().second;
			double cur_distance = -queue.top().first;
			queue.pop();

			if (u == end_vertex)
				return make_pair(build_shortest_path(parent, start_vertex, end_vertex), cur_distance);

			if (cur_distance > distance[u])  continue;

			vector<pair<int, double>> v = graph.get_adjacent_edges(u);
			for (size_t i = 0; i < v.size(); i++)
			{
				if (distance[v[i].first] > v[i].second + distance[u]) {
					distance[v[i].first] = v[i].second + distance[u];
					parent[v[i].first] = u;
					queue.push(make_pair(-distance[v[i].first], v[i].first));
				}
			}
		}
	}

	return make_pair(vector<int> {}, numeric_limits<int>::max());
}

Graph build_visibility_graph(vector<MyPoint> points, size_t rect_size)
{
	Graph graph;
	for (int i = 0; i < points.size(); i++)
		graph.add_vertex(i);

	for (int i = 0; i < points.size() - 1; i++)
		for (int j = i + 1; j < points.size(); j++)
		{
			bool visible = true;
			for (int r = 0; r < rect_size; r++)
			{
				int u_pos = 0;
				int v_pos = 0;
				if (points[i].get_x() <= points[r * 4 + 2].get_x())	u_pos += 1;
				if (points[i].get_x() >= points[r * 4 + 4].get_x())	u_pos += 1 << 1;
				if (points[i].get_y() >= points[r * 4 + 4].get_y())	u_pos += 1 << 2;
				if (points[i].get_y() <= points[r * 4 + 2].get_y())	u_pos += 1 << 3;
				if (points[j].get_x() <= points[r * 4 + 2].get_x())	v_pos += 1;
				if (points[j].get_x() >= points[r * 4 + 4].get_x())	v_pos += 1 << 1;
				if (points[j].get_y() >= points[r * 4 + 4].get_y())	v_pos += 1 << 2;
				if (points[j].get_y() <= points[r * 4 + 2].get_y())	v_pos += 1 << 3;
				if (u_pos & v_pos)
					continue;

				int left = 0, right = 0;
				MySegment seg(points[i], points[j]);
				if (seg.cross(MySegment(points[i], points[r * 4 + 2])) > 0)	left++;
				if (seg.cross(MySegment(points[i], points[r * 4 + 2])) < 0)	right++;
				if (seg.cross(MySegment(points[i], points[r * 4 + 3])) > 0) left++;
				if (seg.cross(MySegment(points[i], points[r * 4 + 3])) < 0) right++;
				if (seg.cross(MySegment(points[i], points[r * 4 + 4])) > 0) left++;
				if (seg.cross(MySegment(points[i], points[r * 4 + 4])) < 0) right++;
				if (seg.cross(MySegment(points[i], points[r * 4 + 5])) > 0) left++;
				if (seg.cross(MySegment(points[i], points[r * 4 + 5])) < 0) right++;
				if (left != 0 && right != 0) {
					visible = false;
					break;
				}
			}
			if (visible)
				graph.add_edge(i, j, points[i].distance(points[j]));
		}

	return graph;
}

pair<vector<int>, float>build_path(Graph& graph, size_t rectangles_size, size_t intermediates_size)
{
	if (intermediates_size)
	{
		vector<int> path;
		double path_lenght = 0;

		vector<int> points_to_access;
		for (int i = int(rectangles_size * 4 + 2); i < 2 + rectangles_size * 4 + intermediates_size; i++)
			points_to_access.push_back(i);

		pair<vector<int>, double> path_to_add = make_pair(vector<int> {}, numeric_limits<int>::max());
		for (auto& e : points_to_access)
		{
			pair<vector<int>, double> new_path = shortest_path(graph, 0, e);
			if (!new_path.first.size())
				return make_pair(vector<int> {}, numeric_limits<int>::max());
			if (new_path.second < path_to_add.second)
				path_to_add = new_path;
		}
		points_to_access.erase(remove(points_to_access.begin(), points_to_access.end(), *(path_to_add.first.end() - 1)), points_to_access.end());
		path.insert(path.end(), path_to_add.first.begin(), path_to_add.first.end());
		path_lenght += path_to_add.second;

		while (points_to_access.size())
		{
			path_to_add = make_pair(vector<int> {}, numeric_limits<int>::max());
			for (auto& e : points_to_access)
			{
				pair<vector<int>, double> new_path = shortest_path(graph, *(path.end() - 1), e);
				if (!new_path.first.size())
					return make_pair(vector<int> {}, numeric_limits<int>::max());
				if (new_path.second < path_to_add.second)
					path_to_add = new_path;
			}
			points_to_access.erase(remove(points_to_access.begin(), points_to_access.end(), *(path_to_add.first.end() - 1)), points_to_access.end());
			path.insert(path.end(), path_to_add.first.begin() + 1, path_to_add.first.end());
			path_lenght += path_to_add.second;
		}

		path_to_add = shortest_path(graph, *(path_to_add.first.end() - 1), 1);
		if (!path_to_add.first.size())
			return path_to_add;
		path.insert(path.end(), path_to_add.first.begin() + 1, path_to_add.first.end());
		path_lenght += path_to_add.second;

		return make_pair(path, path_lenght);
	}
	else
		return shortest_path(graph, 0, 1);
}

pair<pair<vector<MySegment>, vector<MySegment>>, float> build_visibility_graph_and_path(vector<MyRectangle>& rectangles, shared_ptr<MyPoint> start, float ROBOT_SIZE, shared_ptr<MyPoint> finish, vector<MyPoint> intermediates = {})
{
	vector<MyPoint> points;
	points.push_back(*start);
	points.push_back(*finish);
	for (auto& r : rectangles)
	{
		vector<MyPoint> r_points = r.getVertices(ROBOT_SIZE);
		points.insert(points.end(), r_points.begin(), r_points.end());
	}
	for (auto& p : intermediates)
		points.push_back(p);

	Graph graph = build_visibility_graph(points, rectangles.size());

	auto [path, path_lenght] = build_path(graph, rectangles.size(), intermediates.size());

	vector<MySegment> result_graph;
	vector<MySegment> result_path;

	if (path.size())
	{
		for (auto& e : graph.get_edges()) {
			bool visible = true;
			for (int i = 0; i < path.size() - 1; i++)
				if (e.first.first == path[i] && e.first.second == path[i + 1] || e.first.first == path[i + 1] && e.first.second == path[i]) {
					visible = false;
					break;
				}
			if (visible)
				result_graph.push_back(MySegment(points[e.first.first], points[e.first.second]));
		}

		for (int p = 0; p < path.size() - 1; p++)
			result_path.push_back(MySegment(points[path[p]], points[path[p + 1]]));
	}
	else
		for (auto& e : graph.get_edges()) 
			result_graph.push_back(MySegment(points[e.first.first], points[e.first.second]));

	return make_pair(make_pair(result_graph, result_path), path_lenght);
}

void save_scene(string file_name, float& ROBOT_SIZE, shared_ptr<MyPoint>& start, shared_ptr<MyPoint>& finish, vector<MyRectangle>& rectangles, vector<MyPoint>& intermediates)
{
	ofstream file(file_name);

	file << "<ROBOT_SIZE>\n";
	if (start)
		file << ROBOT_SIZE;
	file << "\n</ROBOT_SIZE>\n";

	file << "<start>\n";
	if (start)
		file << *start;
	file << "\n</start>\n";

	file << "<finish>\n";
	if (finish)
		file << *finish;
	file << "\n</finish>\n";

	file << "<rectangles>\n";
	for (auto& r : rectangles)
		file << r << ' ';
	file << "\n</rectangles>\n";

	file << "<intermediates>\n";
	for (auto& p : intermediates)
		file << p << ' ';
	file << "\n</intermediates>\n";
}

bool load_scene(string file_name, float& ROBOT_SIZE, shared_ptr<MyPoint>& start, shared_ptr<MyPoint>& finish, vector<MyRectangle>& rectangles, vector<MyPoint>& intermediates)
{
	ifstream file(file_name);

	float NEW_ROBOT_SIZE = 1.f;
	shared_ptr<MyPoint> new_start;
	shared_ptr<MyPoint> new_finish;
	vector<MyRectangle> new_rectangles;
	vector<MyPoint> new_intermediates;

	if (file.is_open())
	{
		string str;
		float x1, y1, x2, y2;

		file >> str;
		if (str != "<ROBOT_SIZE>")
			return false;

		file >> str;
		if (str != "</ROBOT_SIZE>")
		{
			NEW_ROBOT_SIZE = stof(str);
			file >> str;
			if (str != "</ROBOT_SIZE>")
				return false;
		}

		file >> str;
		if (str != "<start>")
			return false;

		file >> str;
		if (str != "</start>")
		{
			x1 = stof(str);
			file >> y1;
			if (!file)
				return false;
			new_start = make_shared<MyPoint>(x1, y1);
			file >> str;
			if (str != "</start>")
				return false;
		}

		file >> str;
		if (str != "<finish>")
			return false;

		file >> str;
		if (str != "</finish>")
		{
			x1 = stof(str);
			file >> y1;
			if (!file)
				return false;
			new_finish = make_shared<MyPoint>(x1, y1);
			file >> str;
			if (str != "</finish>")
				return false;
		}

		file >> str;
		if (str != "<rectangles>")
			return false;
		
		file >> str;
		while (str != "</rectangles>")
		{
			x1 = stof(str);
			file >> y1;
			file >> x2;
			file >> y2;
			if (!file)
				return false;
			new_rectangles.push_back(MyRectangle(MyPoint(x1, y1), MyPoint(x2, y2)));
			file >> str;
		}

		file >> str;
		if (str != "<intermediates>")
			return false;

		file >> str;
		while (str != "</intermediates>")
		{
			x1 = stof(str);
			file >> y1;
			if (!file)
				return false;
			new_intermediates.push_back(MyPoint(x1, y1));
			file >> str;
		}

		ROBOT_SIZE = NEW_ROBOT_SIZE;
		start = new_start;
		finish = new_finish;
		rectangles = new_rectangles;
		intermediates = new_intermediates;
		return true;
	}

	return false;
}

int main()
{
	float ROBOT_SIZE;
	shared_ptr<MyPoint> start;
	shared_ptr<MyPoint> finish;
	vector<MyRectangle> rectangles;
	vector<MyPoint> intermediates;

	RenderWindow window(VideoMode(1280, 720), "Visibility graph, rectangles, with intermediates.", sf::Style::Titlebar | sf::Style::Close);
	window.setVerticalSyncEnabled(true);
	window.setKeyRepeatEnabled(false);
	Image icon;
	if (!icon.loadFromFile("resources/icon.png"))
		throw invalid_argument("Can't find icon 'icon.png'\n");
	window.setIcon(256, 256, icon.getPixelsPtr());

	Scene scene;

	shared_ptr<MyPoint> firstRectVerticy;

	bool shiftPressed = false;
	bool mouseMoved = false;
	bool graph_linesVisible = false;
	bool buildingMode = false;
	bool helpVisible = true;

	while (window.isOpen())
	{
		Event event;
		while (window.pollEvent(event))
		{
			switch (event.type)
			{
			case Event::Closed:

				window.close();
				break;

			case Event::MouseButtonPressed:

				if (shiftPressed && !buildingMode && event.mouseButton.button == Mouse::Left)
				{
					if (!start)
					{
						start = make_shared<MyPoint>(MyPoint(event.mouseButton.x, event.mouseButton.y));
						ROBOT_SIZE = 5.0f;
					}
					else {
						start->set_x(event.mouseButton.x);
						start->set_y(event.mouseButton.y);
					}
					scene.set_robot(start, ROBOT_SIZE);
					if (finish)
						scene.render_path(build_visibility_graph_and_path(rectangles, start, ROBOT_SIZE, finish, intermediates));
					break;
				}

				if (shiftPressed && !buildingMode && event.mouseButton.button == Mouse::Right)
				{
					if (!finish)
						finish = make_shared<MyPoint>(MyPoint(event.mouseButton.x, event.mouseButton.y));
					else {
						finish->set_x(event.mouseButton.x);
						finish->set_y(event.mouseButton.y);
					};
					scene.set_destination(finish);
					if (start)
						scene.render_path(build_visibility_graph_and_path(rectangles, start, ROBOT_SIZE, finish, intermediates));
					break;
				}

				if (buildingMode && event.mouseButton.button == Mouse::Left)
				{
					if (firstRectVerticy == NULL)
						firstRectVerticy = make_unique<MyPoint>(MyPoint(event.mouseButton.x, event.mouseButton.y));
					else
					{
						MyRectangle new_rect(*firstRectVerticy, MyPoint(event.mouseButton.x, event.mouseButton.y));
						rectangles.push_back(new_rect);
						scene.add_rectangle(new_rect);
						if (start && finish)
							scene.render_path(build_visibility_graph_and_path(rectangles, start, ROBOT_SIZE, finish, intermediates));
						scene.upd_preview_rect();
						firstRectVerticy = NULL;
						mouseMoved = false;
					}
					break;
				}

				if (buildingMode && !firstRectVerticy && event.mouseButton.button == Mouse::Right && rectangles.size())
				{
					MyPoint place(event.mouseButton.x, event.mouseButton.y);
					MyRectangle r_to_delete = *place.find_closest_rectangle(rectangles);
					if (place.in(r_to_delete))
					{
						scene.remove_rectangle(r_to_delete);
						rectangles.erase(remove(rectangles.begin(), rectangles.end(), r_to_delete), rectangles.end());
						if (start && finish)
							scene.render_path(build_visibility_graph_and_path(rectangles, start, ROBOT_SIZE, finish, intermediates));
					}
					break;
				}

				if (firstRectVerticy && event.mouseButton.button == Mouse::Right)
				{
					firstRectVerticy = NULL;
					mouseMoved = false;
					scene.upd_preview_rect();
					break;
				}

				if (event.mouseButton.button == Mouse::Left)
				{
					MyPoint p_to_add(event.mouseButton.x, event.mouseButton.y);
					if (!intermediates.size() || p_to_add.distance(*p_to_add.find_closest(intermediates)) > 0)
					{
						intermediates.push_back(p_to_add);
						scene.add_dot(intermediates.back());
						if (start && finish)
							scene.render_path(build_visibility_graph_and_path(rectangles, start, ROBOT_SIZE, finish, intermediates));
					}
				}

				if (event.mouseButton.button == Mouse::Right && intermediates.size())
				{
					MyPoint place(event.mouseButton.x, event.mouseButton.y);
					MyPoint p_to_delete = *place.find_closest(intermediates);
					if (place.distance(p_to_delete) <= 5)
					{
						scene.remove_dot(p_to_delete);
						intermediates.erase(remove(intermediates.begin(), intermediates.end(), p_to_delete), intermediates.end());
						if (start && finish)
							scene.render_path(build_visibility_graph_and_path(rectangles, start, ROBOT_SIZE, finish, intermediates));
					}
				}
				break;

			case Event::MouseButtonReleased:

				if (firstRectVerticy && mouseMoved && event.mouseButton.button == Mouse::Left)
				{
					MyRectangle new_rect(*firstRectVerticy, MyPoint(event.mouseButton.x, event.mouseButton.y));
					rectangles.push_back(new_rect);
					scene.add_rectangle(new_rect);
					if (start && finish)
						scene.render_path(build_visibility_graph_and_path(rectangles, start, ROBOT_SIZE, finish, intermediates));
					scene.upd_preview_rect();
					firstRectVerticy = NULL;
					mouseMoved = false;
				}

				break;

			case Event::MouseMoved:

				if (firstRectVerticy != NULL)
				{
					mouseMoved = true;
					scene.upd_preview_rect(firstRectVerticy, event.mouseMove.x, event.mouseMove.y);
				}

				break;

			case Event::MouseWheelScrolled:

				if (shiftPressed && !buildingMode && event.mouseWheelScroll.wheel == Mouse::VerticalWheel)
				{
					ROBOT_SIZE += ROBOT_SIZE + event.mouseWheelScroll.delta > 0 ? event.mouseWheelScroll.delta : -ROBOT_SIZE;
					if (start)
						scene.set_robot(start, ROBOT_SIZE);
					if (start && finish)
						scene.render_path(build_visibility_graph_and_path(rectangles, start, ROBOT_SIZE, finish, intermediates));
				}

				break;

			case Event::KeyReleased:

				if (event.key.code == sf::Keyboard::LShift || event.key.code == sf::Keyboard::RShift)
					shiftPressed = false;

				break;

			case Event::KeyPressed:

				if (buildingMode && event.key.code == sf::Keyboard::Escape)
				{
					firstRectVerticy = NULL;
					mouseMoved = false;
					buildingMode = false;
					scene.upd_preview_rect();
					break;
				}

				if (event.key.code == sf::Keyboard::Escape)
				{
					window.close();
					break;
				}

				if (event.key.code == sf::Keyboard::F5)
				{
					save_scene("saves/quicksave", ROBOT_SIZE, start, finish, rectangles, intermediates);
					break;
				}

				if (event.key.code == sf::Keyboard::F9)
				{
					if (load_scene("saves/quicksave", ROBOT_SIZE, start, finish, rectangles, intermediates))
					{

						scene.clear();
						if (start) scene.set_robot(start, ROBOT_SIZE);
						if (finish) scene.set_destination(finish);
						for (auto& r : rectangles) scene.add_rectangle(r);
						for (auto& d : intermediates) scene.add_dot(d);
						if (start && finish)
							scene.render_path(build_visibility_graph_and_path(rectangles, start, ROBOT_SIZE, finish, intermediates));
					}
					break;
				}

				if (event.key.code == sf::Keyboard::LShift || event.key.code == sf::Keyboard::RShift)
					shiftPressed = true;

				if (event.key.code == sf::Keyboard::H)
				{
					if (graph_linesVisible)
						graph_linesVisible = false;
					else
						graph_linesVisible = true;
				}

				if (event.key.code == sf::Keyboard::B)
				{
					if (buildingMode)
					{
						firstRectVerticy = NULL;
						mouseMoved = false;
						buildingMode = false;
						scene.upd_preview_rect();
					}
					else
						buildingMode = true;
				}

				if (event.key.code == sf::Keyboard::F1)
				{
					if (helpVisible)
						helpVisible = false;
					else
						helpVisible = true;
				}

				if (shiftPressed && buildingMode && event.key.code == sf::Keyboard::C)
				{
					rectangles.clear();
					scene.rectangles.clear();
					if (start && finish)
						scene.render_path(build_visibility_graph_and_path(rectangles, start, ROBOT_SIZE, finish, intermediates));
					break;
				}

				if (shiftPressed && event.key.code == sf::Keyboard::C)
				{
					intermediates.clear();
					scene.dots.clear();
					if (start && finish)
						scene.render_path(build_visibility_graph_and_path(rectangles, start, ROBOT_SIZE, finish, intermediates));
				}

				if (shiftPressed && event.key.code == sf::Keyboard::D)
				{
					start = NULL;
					finish = NULL;
					intermediates.clear();
					rectangles.clear();
					scene.clear();
				}

				break;

			default:
				break;
			}
		}

		window.clear(sf::Color(240, 248, 255));

		for (auto& rect : scene.rectangles)
			window.draw(rect);

		if (firstRectVerticy != NULL)
			window.draw(scene.preview_rect);

		if (graph_linesVisible)
			for (auto& line : scene.graph_lines)
				window.draw(line);

		for (auto& line : scene.path_lines)
			window.draw(line);

		for (auto& dot : scene.dots)
			window.draw(dot);

		window.draw(scene.robot);
		window.draw(scene.destination);

		window.draw(scene.path_lenght);

		if (buildingMode)
			window.draw(scene.edit_mode_indicator);

		if (helpVisible)
		{
			window.draw(scene.help_underlay);
			window.draw(scene.help);
		}

		window.display();
	}

	return EXIT_SUCCESS;
}