#ifndef CONSTANTS_HPP
# define CONSTANTS_HPP

enum e_tags {
  NONE       = 0,
  FOREGROUND = 1,  // Object
  BACKGROUND = 2,  // Background
  DUMMY      = 3,  // 3rd class
  REJECTED   = 4  // Rejection label
};

const unsigned int AMBIGUITY_MASK = 0x07;

enum e_labelization_policy {
  AMBIGUITY_AS_FOREGROUND = 0x00,
  AMBIGUITY_AS_BACKGROUND = 0x01,
  AMBIGUITY_AS_DUMMY      = 0x02,
  AMBIGUITY_AS_MAJORITY   = 0x03,
  AMBIGUITY_AS_UNTAGGED   = 0x04
};


#endif // ! CONSTANTS_HPP
