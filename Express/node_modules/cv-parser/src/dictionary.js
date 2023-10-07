const dictionary = {
  titles: {
    objective: ["objective", "objectives"],
    summary: ["summary"],
    technology: ["technology", "technologies"],
    experience: ["experience", "experiences"],
    volunteering: ["volunteering", "volunteer", "probono"],
    education: ["education"],
    skills: [
      "skills",
      "Skills & Expertise",
      "technologies",
      "technical skills",
    ],
    languages: ["languages"],
    courses: ["courses"],
    projects: ["projects"],
    links: ["links"],
    contacts: ["contacts"],
    positions: ["positions", "position"],
    profiles: [
      "profiles",
      "social connect",
      "social-profiles",
      "social profiles",
    ],
    awards: ["awards"],
    honors: ["honors"],
    additional: ["additional"],
    certification: ["certification", "certifications"],
    interests: ["interests"],
  },
  regular: {
    name: [/([A-Z][a-z]*)(\s[A-Z][a-z]*)/],
    email: [/([a-z0-9_\.-]+)@([\da-z\.-]+)\.([a-z\.]{2,6})/],
    phone: [/^(\+\d{1,5}\s)?\(?\d{3}\)?[\s.-]\d{3}[\s.-]\d{4}$/],
  },
  profiles: ["github.com"],
};

module.exports = dictionary;
