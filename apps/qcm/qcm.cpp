#include <mln/core/image/image2d.hpp>
#include <mln/core/colors.hpp>
#include <mln/core/algorithm/transform.hpp>
#include <mln/core/algorithm/accumulate.hpp>
#include <mln/core/algorithm/copy.hpp>
#include <mln/accu/accumulators/sum.hpp>
#include <mln/io/imread.hpp>
#include <mln/io/imsave.hpp>

#include <mln/core/win2d.hpp>
#include <mln/morpho/structural/opening.hpp>
#include <mln/core/image/morphers/casted_image.hpp>
#include <boost/program_options.hpp>


namespace mln
{

  static const     rect2d qcmbox = make_rectangle2d(15, 41);
  static constexpr float qcmboxvdist = 50.2;
  static constexpr int qcmboxhdist = 60;
  static constexpr int qcmcolhdist = 481;
  static constexpr point2d pref = {168, 2389};

  static constexpr char REPA = 0x01;
  static constexpr char REPB = 0x02;
  static constexpr char REPC = 0x04;
  static constexpr char REPD = 0x08;
  static constexpr char REPE = 0x10;
  static constexpr char REPALL = 0x1F;
  static unsigned detection_threshold = 200;


  image2d<rgb8> imdebug;




  image2d<bool>
  transpose(const image2d<bool>& x)
  {
    box2d b1 = x.domain();
    box2d b2 = {{b1.pmin[1], b1.pmin[0]},{b1.pmax[1], b1.pmax[0]}};

    image2d<bool> out(b2);


    mln_foreach(auto pxin, x.pixels())
    {
      point2d p = pxin.point();
      out.at(p[1],p[0]) = pxin.val();
    }

    return out;
  }

  point2d
  detect_offset(const image2d<uint8>& f)
  {
    mln_entering("Offset detection");
    //rect2d se = make_rectangle2d(21, 51);

    box2d dom = { {125, 2300}, {250, 2470} };
    image2d<bool> bin(dom);
    copy((f < 150) | dom, bin);

    //auto markers = morpho::opening(f, se);
    rect2d se = make_rectangle2d(21,51);
    bin = morpho::structural::opening(bin, se);

    //auto sub = bin | dom;
    mln_foreach(auto px, bin.pixels())
      if (px.val()) {
        mln_exiting();
        return px.point();
      }

    std::cerr << "Marker not found" << std::endl;
    std::exit(1);
  }

  template <class I>
  bool
  is_plain(const I& f, point2d p)
  {
    int count = accumulate(f | qcmbox(p), accu::features::sum<>());

    if (count > detection_threshold)
      {
        for (int i = p[1] - 20; i <= (p[1] + 20); ++i)
          {
            imdebug.at(p[0]-8,i) = rgb8{255,0,0};
            imdebug.at(p[0]-7,i) = rgb8{255,0,0};
            imdebug.at(p[0]+7,i) = rgb8{255,0,0};
            imdebug.at(p[0]+8,i) = rgb8{255,0,0};
          }
      }

    //std::cout << p << ":" << count << std::endl;
    return count > detection_threshold;
  }

  template <class I>
  char
  detect_question(const I& bin, point2d pos)
  {
    char response1 = 0;
    char response2 = 0;
    point2d p = pos;

    for (int i = 0; i < 5; ++i, p[1] += qcmboxhdist)
      if (is_plain(bin, p))
        response1 |= (1 << i);

    p = pos;
    p[0] += qcmboxvdist;

    for (int i = 0; i < 5; ++i, p[1] += qcmboxhdist)
      if (is_plain(bin, p))
        response2 |= (1 << i);

    if (response2 == REPALL)
      return 0;
    else if (response2)
      return response2;
    else if (response1 and response1 != REPALL)
      return response1;
    else
      return 0;
  }

  void
  detect_all(const image2d<uint8>& f, point2d offset = {0,0})
  {
    mln_entering("Detection des cases");

    auto bin = f < 200;

    rect2d se = make_rectangle2d(15, 40);

    // Login
    {
      point2d plogin = point2d{228, 1895} + offset;
      point2df p = plogin;
      for (int j = 0; j < 6; ++j, p[1] += qcmboxhdist)
        {
          p[0] = plogin[0];
          for (int i = 0; i < 26; ++i, p[0] += qcmboxvdist)
            {
              if (is_plain(bin, p)) {
                std::cout << (char)('a' + i);
                break;
              }
            }
        }
      std::cout << "-";
      p[1] += qcmboxhdist;
      p[0] = plogin[0];
      for (int i = 0; i < 26; ++i, p[0] += qcmboxvdist)
        {
          if (is_plain(bin, p)) {
            std::cout << (char)('a' + i);
            break;
          }
        }
      std::cout << std::endl;
    }


    // Les questions
    {
      point2d p0 = point2d{1734, 156} + offset;

      point2df p = p0;
      for (int j = 0; j < 5; ++j, p[1] += qcmcolhdist) {
        p[0] = p0[0];
        for (int i = 0; i < 10; ++i, p[0] += 2*qcmboxvdist) {
          char res = detect_question(bin, p);
          //std::cout << (j*10 + i + 1) << " ";

          for (int k = 0; k < 5; ++k)
            if (res & (1 << k))
              std::cout << (char)('a' + k);

          std::cout << std::endl;
        }
      }
    }
    mln_exiting();
  }


}




int main(int argc, char** argv)
{
  using namespace mln;
  namespace po = boost::program_options;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("threshold,t", po::value<unsigned>(&detection_threshold)->default_value(200), "Set threshold detection")
    ;

  po::options_description desc2("");
  desc2.add_options()
    ("input", po::value<std::string>()->required(), "input image")
    ("output", po::value<std::string>()->required(), "output image");

  po::positional_options_description pd;
  pd.add("input", 1);
  pd.add("output", 1);

  desc2.add(desc);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(desc2).positional(pd).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  image2d<uint8> f;

  {
    mln_entering("Conversion B&W")
    image2d<rgb8> ima;
    io::imread(vm["input"].as<std::string>(), ima);

    f = transform(ima, [](const rgb8& v) -> uint8 { return sum(v) / 3;});
    imdebug = ima;
    copy(f, imdebug);
    mln_exiting();
  }

  point2d ref = detect_offset(f);
  point2d offset = ref - pref;

  detect_all(f, offset);
  io::imsave(imdebug, vm["output"].as<std::string>());

  // std::cout << "Ref: " << ref << std::endl;
  // std::cout << "Offset: " << offset << std::endl;
}
