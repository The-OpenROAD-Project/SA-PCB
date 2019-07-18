class Pin {
  public:
  string name;
  double x_offset;
  double y_offset;
  int idx;
  double depth;

  void set_params(string nn, double xoffset, double yoffset, int idx) {
    this -> name = nn;
    this -> x_offset = xoffset;
    this -> y_offset = yoffset;
    this -> idx = idx;
  }
};
