---
# Leave the homepage title empty to use the site title
title: ''
summary: ''
date: 2026-04-13
type: landing

design:
  # Default section spacing
  spacing: '6rem'

sections:
  - block: resume-biography-3
    content:
      # Choose a user profile to display (a folder name within `content/authors/`)
      username: me
      text: ''
      headings:
        about: ''
        education: ''
        interests: ''
    design:
      # Use the new Gradient Mesh which automatically adapts to the selected theme colors
      background:
        gradient_mesh:
          enable: true

      # Name heading sizing to accommodate long or short names
      name:
        size: md # Options: xs, sm, md, lg (default), xl

      # Avatar customization
      avatar:
        size: medium # Options: small (150px), medium (200px, default), large (320px), xl (400px), xxl (500px)
        shape: circle # Options: circle (default), square, rounded
  - block: markdown
    content:
      title: 'About Me'
      subtitle: ''
      text: |-
        I'm a nerd for computers and biology. I've tried to bring those two things together throughout my PhD. I like to
        tackle problems from multiple angles, and like learning how to work with new data. I've worked with a number
        of aspects of scRNA-seq data, from mutations, to time inference, to analyzing RNA counts. In my spare time I like
        reading papers about disease management and health outcomes. Directly applying my work to this would be quite
        exciting for me, so please contact me if you have any opportunities!
        
        Aside from my work, I like to go outside ... a lot. I like hiking, mountaineering, climbing, and everything in-between.
        Lately, I've also found a bit of a niche in maintaining backcountry huts. I recently organized two large hut
        restorations with the UBC Varsity Outdoor Club, which brought together hundreds of volunteers over thousands of
        hours of work to save two aging huts. I post trip reports about escapades like this and more on this site.
    design:
      columns: '1'  
  - block: collection
    content:
      title: Recent Publications
      text: ''
      filters:
        folders:
          - publications
        exclude_featured: false
    design:
      view: citation
---
